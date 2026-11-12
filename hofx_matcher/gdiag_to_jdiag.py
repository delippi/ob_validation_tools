#!/usr/bin/env python3
import argparse
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
import os

# Platform mapping inferred from jdiag file list
platform_map = {
    120: "adpupa", 132: "adpupa", 220: "adpupa", 232: "adpupa",
    130: "aircft", 131: "aircft", 134: "aircft", 135: "aircft",
    230: "aircft", 231: "aircft", 234: "aircft", 235: "aircft",
    133: "aircar", 233: "aircar",
    180: "sfcshp", 182: "sfcshp", 183: "sfcshp",
    280: "sfcshp", 282: "sfcshp", 284: "sfcshp",
    181: "adpsfc", 183: "adpsfc", 187: "adpsfc",
    281: "adpsfc", 284: "adpsfc", 287: "adpsfc",
    188: "msonet", 288: "msonet",
    224: "vadwnd",
    227: "proflr",
    126: "rassda",
    242: "satwnd",
    245: "satwnd", 246: "satwnd", 247: "satwnd",
    250: "satwnd", 252: "satwnd", 253: "satwnd",
    999: "nexrad",
}

# Variable mapping from GSI to JEDI
variable_map = {
    "t": "airTemperature",
    "q": "specificHumidity",
    "ps": "stationPressure",
    "uv": "winds",  # Special case for windEastward and windNorthward
    "rw": "radialVelocity",
}

def gsi_time_to_epoch(gsi_time, analysis_time_str):
    """Convert GSI time offset to epoch time in seconds since 1970-01-01."""
    analysis_time = datetime.strptime(analysis_time_str, '%Y%m%d%H')
    obs_time = analysis_time + timedelta(hours=float(gsi_time))
    epoch = datetime(1970, 1, 1)
    return int((obs_time - epoch).total_seconds())

def char_array_to_strings(char_array):
    """
    Convert a 2D character array to a list of strings, handling masked arrays.

    Args:
        char_array: 2D numpy array or masked array of characters (nobs, maxstrlen)

    Returns:
        List of strings, with masked rows replaced by 'MISSING'
    """
    strings = []
    for row in char_array:
        # GSI radar IDs contain four characters in a longer character field.
        # NetCDF4 returns the unused trailing positions as MaskedConstant values.
        masked_row = np.ma.asarray(row)
        values = masked_row.compressed() if masked_row.ndim else [masked_row.item()]

        characters = []
        for value in values:
            if np.ma.is_masked(value):
                continue
            if isinstance(value, (bytes, np.bytes_)):
                characters.append(value.decode('ascii', errors='ignore'))
            else:
                characters.append(str(value))

        # Remove trailing blanks and null characters from the fixed-width field.
        station_id = ''.join(characters).replace('\x00', '').strip()
        strings.append(station_id if station_id else 'MISSING')
    return strings

def inverse_observation_error(errinv):
    """Convert inverse errors to errors without dividing missing or zero values."""
    errinv_array = np.ma.asarray(errinv)
    errinv_data = np.asarray(np.ma.getdata(errinv_array), dtype=np.float64)
    errinv_mask = np.ma.getmaskarray(errinv_array)

    # Values below the smallest normal float32 are not usable inverse errors.
    valid = (
        ~errinv_mask
        & np.isfinite(errinv_data)
        & (errinv_data > np.finfo(np.float32).tiny)
    )

    effective_error = np.ma.masked_all(errinv_data.shape, dtype=np.float32)
    effective_error[valid] = (1.0 / errinv_data[valid]).astype(np.float32)
    return effective_error

def map_prep_use_flag_to_qc(flag):
    """Map GSI Prep_Use_Flag to JEDI EffectiveQC values."""
    if flag == 1.0:
        return 0  # Assimilated
    elif flag == -1.0:
        return 1  # Monitored
    else:
        return 12  # Rejected

def process_non_wind_file(
    anl_file_path, ges_file_path, base_jedi_var, analysis_time_str, surface_at2m=False
):
    """
    Process a non-wind GSI diagnostic file and convert it to JEDI jdiag format.

    Args:
        file_path (str): Path to the GSI diagnostic netCDF file
        base_jedi_var (str): JEDI variable name before platform-specific naming
        analysis_time_str (str): Analysis time string for time conversion
        surface_at2m (bool): Use At2M names for surface temperature and humidity
    """
    with nc.Dataset(ges_file_path, 'r') as ges_f, nc.Dataset(anl_file_path, 'r') as anl_f:
        # Read variables from the GSI file
        station_id_char = ges_f.variables['Station_ID'][:]  # 2D char array (nobs, maxstrlen)
        obs_type = ges_f.variables['Observation_Type'][:]
        latitude = ges_f.variables['Latitude'][:]
        longitude = ges_f.variables['Longitude'][:]
        pressure = ges_f.variables['Pressure'][:]
        # Radar radial-velocity diagnostics use geometric height in JEDI.
        height = ges_f.variables['Height'][:] if 'Height' in ges_f.variables else None
        station_elevation = (
            ges_f.variables['Station_Elevation'][:]
            if 'Station_Elevation' in ges_f.variables else None
        )
        time = ges_f.variables['Time'][:]
        #effective_error_0 = np.sqrt(1.0 / ges_f.variables['Errinv_Final'][:])  # Invert
        effective_error_0 = inverse_observation_error(ges_f.variables['Errinv_Final'][:])
        observation = ges_f.variables['Observation'][:]
        ombg = ges_f.variables['Obs_Minus_Forecast_adjusted'][:]
        oman = anl_f.variables['Obs_Minus_Forecast_adjusted'][:]
        analysis_use_flag = anl_f.variables['Analysis_Use_Flag'][:]
        prep_use_flag = anl_f.variables['Prep_Use_Flag'][:]
        # Convert stationPressure units: GSI uses hPa, JEDI expects Pa
        if base_jedi_var == "stationPressure":
            observation = observation * 100.0
            ombg = ombg * 100.0
            oman = oman * 100.0
            effective_error_0 = effective_error_0 * 100.0

        # Process each unique observation type
        unique_types = np.unique(obs_type)
        for typ in unique_types:
            if typ not in platform_map:
                print(f"Warning: Observation_Type {typ} not in platform_map, skipping.")
                continue
            platform = platform_map[typ]
            # One GSI file can contain both surface and upper-air platforms.
            jedi_var = base_jedi_var
            if (
                surface_at2m
                and platform in ("adpsfc", "msonet", "sfcshp")
                and base_jedi_var in ("airTemperature", "specificHumidity")
            ):
                jedi_var += "At2M"
            mask = obs_type == typ
            nobs = int(np.sum(mask))
            if nobs == 0:
                print(f"No observations for type {typ} in file {anl_file_path}, skipping.")
                continue

            # Convert selected station IDs to strings
            selected_station_ids = char_array_to_strings(station_id_char[mask])

            # Use combination of Prep_Use and Analysis_Use flags to determine asm, rej, mon
            # Preserve one QC value per observation even if a flag is masked.
            iuse_valid = np.ma.filled(prep_use_flag[mask], 0.0)
            ause_valid = np.ma.filled(analysis_use_flag[mask], 0.0)

            # Initialize flag array with zeros (default for rejected)
            flag = np.zeros_like(ause_valid, dtype=float)

            # Set flag for assimilated observations
            flag[ause_valid == 1] = 1.0

            # Set flag for monitored observations
            flag[(ause_valid == -1) & (iuse_valid > 0)] = -1.0

            # Map to JEDI EffectiveQC2 equivalent
            qc_flags = np.array([map_prep_use_flag_to_qc(f) for f in flag], dtype=np.int32)

            # Create output jdiag file
            jdiag_file = f"jdiag_{platform}_{jedi_var}_{typ}.nc"
            with nc.Dataset(jdiag_file, 'w', format='NETCDF4') as g:
                g.setncattr('_ioda_layout', 'ObsGroup')
                g.setncattr('_ioda_layout_version', 0)
                g.createDimension('Location', nobs)

                # MetaData group
                meta = g.createGroup('MetaData')
                # Create the variable for station IDs
                station_var = meta.createVariable('stationIdentification', 'str', ('Location',))
                # Assign each string individually to avoid VLEN assignment error
                for i, s in enumerate(selected_station_ids):
                    station_var[i] = s
                meta.createVariable('latitude', 'f4', ('Location',))[:] = latitude[mask]
                meta.createVariable('longitude', 'f4', ('Location',))[:] = longitude[mask]
                if jedi_var == "radialVelocity":
                    if height is None:
                        raise KeyError(
                            f"Required GSI variable 'Height' is missing from {ges_file_path}"
                        )
                    # Match the vertical coordinate in the NEXRAD JEDI diagnostic.
                    meta.createVariable('height', 'f4', ('Location',))[:] = height[mask]
                    if station_elevation is not None:
                        meta.createVariable(
                            'stationElevation', 'f4', ('Location',)
                        )[:] = station_elevation[mask]
                else:
                    meta.createVariable('pressure', 'f4', ('Location',))[:] = pressure[mask] * 100  # hPa to Pa
                #date_times = [gsi_time_to_epoch(t, analysis_time_str) for t in time[mask]]
                #meta.createVariable('dateTime', 'i8', ('Location',))[:] = date_times
                # dateTime (epoch seconds)
                date_times = [gsi_time_to_epoch(t, analysis_time_str) for t in time[mask]]
                meta.createVariable('dateTime', 'i8', ('Location',))[:] = date_times

                # timeOffset (GSI "Time" is hours relative to analysis)
                time_offset = time[mask] * 3600.0
                meta.createVariable('timeOffset', 'f4', ('Location',))[:] = time_offset


                # ObsValue group
                obs_val = g.createGroup('ObsValue')
                obs_val.createVariable(jedi_var, 'f4', ('Location',))[:] = observation[mask]

                # Retain the GSI observation type in the radial-velocity output.
                if jedi_var == "radialVelocity":
                    obs_type_group = g.createGroup('ObsType')
                    obs_type_group.createVariable(
                        jedi_var, 'i4', ('Location',)
                    )[:] = np.asarray(obs_type[mask], dtype=np.int32)

                # EffectiveError0 group
                eff_err_0 = g.createGroup('EffectiveError0')
                eff_err_0.createVariable(jedi_var, 'f4', ('Location',))[:] = effective_error_0[mask]

                # hofx0 (background)
                hofx0 = g.createGroup('hofx0')
                hofx0.createVariable(jedi_var, 'f4', ('Location',))[:] = observation[mask] - ombg[mask]

                # hofx1 (analysis)
                hofx1 = g.createGroup('hofx1')
                hofx1.createVariable(jedi_var, 'f4', ('Location',))[:] = observation[mask] - oman[mask]

                # ombg group
                ombg_group = g.createGroup('ombg')
                ombg_group.createVariable(jedi_var, 'f4', ('Location',))[:] = ombg[mask]

                # oman group
                oman_group = g.createGroup('oman')
                oman_group.createVariable(jedi_var, 'f4', ('Location',))[:] = oman[mask]

                # EffectiveQC2 group
                eff_qc2 = g.createGroup('EffectiveQC2')
                eff_qc2.createVariable(jedi_var, 'i4', ('Location',))[:] = qc_flags

def process_wind_file(anl_file_path, ges_file_path, analysis_time_str):
    """
    Process a wind GSI diagnostic file and convert it to JEDI jdiag format.

    Args:
        file_path (str): Path to the GSI diagnostic netCDF file
        analysis_time_str (str): Analysis time string for time conversion

    Notes:
        Assumes platform_map and gsi_time_to_epoch are defined elsewhere.
    """
    with nc.Dataset(ges_file_path, 'r') as ges_f, nc.Dataset(anl_file_path, 'r') as anl_f:
        # Read variables from the GSI file
        station_id_char = ges_f.variables['Station_ID'][:]  # 2D char array (nobs, maxstrlen)
        obs_type = ges_f.variables['Observation_Type'][:]
        obs_subtype = ges_f.variables["Observation_Subtype"][:]
        latitude = ges_f.variables['Latitude'][:]
        longitude = ges_f.variables['Longitude'][:]
        pressure = ges_f.variables['Pressure'][:]
        time = ges_f.variables['Time'][:]
        effective_error_0 = inverse_observation_error(ges_f.variables['Errinv_Final'][:])
        u_observation = ges_f.variables['u_Observation'][:]
        v_observation = ges_f.variables['v_Observation'][:]
        u_ombg = ges_f.variables['u_Obs_Minus_Forecast_adjusted'][:]
        v_ombg = ges_f.variables['v_Obs_Minus_Forecast_adjusted'][:]
        u_oman = anl_f.variables['u_Obs_Minus_Forecast_adjusted'][:]
        v_oman = anl_f.variables['v_Obs_Minus_Forecast_adjusted'][:]
        analysis_use_flag = anl_f.variables['Analysis_Use_Flag'][:]
        prep_use_flag = anl_f.variables['Prep_Use_Flag'][:]

        # Process each unique observation type
        unique_types = np.unique(obs_type)
        for typ in unique_types:
            if typ not in platform_map:
                print(f"Warning: Observation_Type {typ} not in platform_map, skipping.")
                continue
            platform = platform_map[typ]
            mask_type = (obs_type == typ)  # Boolean mask for this obs_type

            unique_subtypes = np.unique(obs_subtype[mask_type])
            for sub in unique_subtypes:
                mask = mask_type & (obs_subtype == sub)
                nobs = int(np.sum(mask))
                if nobs == 0:
                    print(f"No obs for type {typ} in file {anl_file_path}, skipping.")
                    continue

                # Select station IDs for u components, limited to nobs
                selected_station_ids = char_array_to_strings(station_id_char[mask])

                # Use combination of Prep_Use and Analysis_Use flags to determine asm, rej, mon
                # Preserve one QC value per observation even if a flag is masked.
                iuse_valid = np.ma.filled(prep_use_flag[mask], 0.0)
                ause_valid = np.ma.filled(analysis_use_flag[mask], 0.0)

                # Initialize flag array with zeros (default for rejected)
                flag = np.zeros_like(ause_valid, dtype=float)

                # Set flag for assimilated observations
                flag[ause_valid == 1] = 1.0

                # Set flag for monitored observations
                flag[(ause_valid == -1) & (iuse_valid > 0)] = -1.0

                # Map to JEDI EffectiveQC2 equivalent
                qc_flags = np.array([map_prep_use_flag_to_qc(f) for f in flag], dtype=np.int32)

                # Create output jdiag file
                if platform == "satwnd":
                    jdiag_file = f"jdiag_{platform}_winds_{typ}_{sub}.nc"
                else:
                    jdiag_file = f"jdiag_{platform}_winds_{typ}.nc"
                with nc.Dataset(jdiag_file, 'w', format='NETCDF4') as g:
                    g.setncattr('_ioda_layout', 'ObsGroup')
                    g.setncattr('_ioda_layout_version', 0)
                    g.createDimension('Location', nobs)

                    # MetaData group
                    meta = g.createGroup('MetaData')
                    # Create the variable for station IDs
                    station_var = meta.createVariable('stationIdentification', 'str', ('Location',))
                    # Assign each string individually to avoid VLEN assignment error
                    for i, s in enumerate(selected_station_ids):
                        station_var[i] = s
                    meta.createVariable('latitude', 'f4', ('Location',))[:] = latitude[mask]
                    meta.createVariable('longitude', 'f4', ('Location',))[:] = longitude[mask]
                    meta.createVariable('pressure', 'f4', ('Location',))[:] = pressure[mask] * 100  # hPa to Pa
                    date_times = [gsi_time_to_epoch(t, analysis_time_str) for t in time[mask]]
                    meta.createVariable('dateTime', 'i8', ('Location',))[:] = date_times
                    time_offset = time[mask] * 3600.0
                    meta.createVariable('timeOffset', 'f4', ('Location',))[:] = time_offset


                    # ObsValue group
                    obs_val = g.createGroup('ObsValue')
                    obs_val.createVariable('windEastward', 'f4', ('Location',))[:] = u_observation[mask]
                    obs_val.createVariable('windNorthward', 'f4', ('Location',))[:] = v_observation[mask]

                    # EffectiveError0 group
                    eff_err_0 = g.createGroup('EffectiveError0')
                    eff_err_0.createVariable('windEastward', 'f4', ('Location',))[:] = effective_error_0[mask]
                    eff_err_0.createVariable('windNorthward', 'f4', ('Location',))[:] = effective_error_0[mask]

                    # hofx0 (background)
                    hofx0 = g.createGroup('hofx0')
                    hofx0.createVariable('windEastward', 'f4', ('Location',))[:] = u_observation[mask] - u_ombg[mask]
                    hofx0.createVariable('windNorthward', 'f4', ('Location',))[:] = v_observation[mask] - v_ombg[mask]

                    # hofx1 (analysis)
                    hofx1 = g.createGroup('hofx1')
                    hofx1.createVariable('windEastward', 'f4', ('Location',))[:] = u_observation[mask] - u_oman[mask]
                    hofx1.createVariable('windNorthward', 'f4', ('Location',))[:] = v_observation[mask] - v_oman[mask]

                    # ombg group
                    ombg_group = g.createGroup('ombg')
                    ombg_group.createVariable('windEastward', 'f4', ('Location',))[:] = u_ombg[mask]
                    ombg_group.createVariable('windNorthward', 'f4', ('Location',))[:] = v_ombg[mask]

                    # oman group
                    oman_group = g.createGroup('oman')
                    oman_group.createVariable('windEastward', 'f4', ('Location',))[:] = u_oman[mask]
                    oman_group.createVariable('windNorthward', 'f4', ('Location',))[:] = v_oman[mask]
                    # EffectiveQC2 group
                    eff_qc2 = g.createGroup('EffectiveQC2')
                    eff_qc2.createVariable('windEastward', 'i4', ('Location',))[:] = qc_flags
                    eff_qc2.createVariable('windNorthward', 'i4', ('Location',))[:] = qc_flags

def convert_gsi_to_jdiag(file_list, analysis_time_str, surface_at2m=False):
    """Convert a list of GSI diagnostic files to JEDI jdiag format."""
    red = "\x1b[31m"
    green = "\x1b[32m"
    normal = "\x1b[0m"
    for file_path in file_list:
        if '_anl' not in file_path:
            #print(f"{red}Skipping non-analysis file: {file_path}{normal}")
            continue
        filename = os.path.basename(file_path)
        parts = filename.split('_')
        if len(parts) < 4 or parts[1] != 'conv':
            #print(f"{red}Skipping invalid file: {filename}{normal}")
            continue
        var = parts[2]  # e.g., 't', 'uv', 'q', 'ps'
        stage = parts[3].split('.')[0]  # Should be 'anl'
        if stage != 'anl':
            #print(f"{red}Skipping non-anl file: {filename}{normal}")
            continue
        if var not in variable_map:
            #print(f"{red}Skipping unhandled variable: {var} in {filename}{normal}")
            continue
        jedi_var = variable_map[var]
        anl_file_path = file_path
        ges_file_path = file_path.replace("_anl", "_ges")
        # Fail with a direct message instead of a later NetCDF open error.
        if not os.path.isfile(anl_file_path):
            raise FileNotFoundError(f"Analysis diagnostic not found: {anl_file_path}")
        if not os.path.isfile(ges_file_path):
            raise FileNotFoundError(f"Guess diagnostic not found: {ges_file_path}")
        if var == 'uv':
            process_wind_file(anl_file_path, ges_file_path, analysis_time_str)
        else:
            process_non_wind_file(
                anl_file_path, ges_file_path, jedi_var, analysis_time_str, surface_at2m
            )
        anl_filename = filename
        ges_filename = filename.replace("_anl", "_ges")
        print(f"{green}Processed {ges_filename} & {anl_filename} into jdiag files.{normal}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GSI diagnostics to JEDI format.")
    parser.add_argument("analysis_time", metavar="YYYYMMDDHH")
    parser.add_argument("files", metavar="ANALYSIS_FILE", nargs="+")
    parser.add_argument(
        "--surface-at2m", action="store_true",
        help="Use At2M names for temperature and humidity from adpsfc, msonet, and sfcshp.",
    )
    args = parser.parse_args()
    convert_gsi_to_jdiag(args.files, args.analysis_time, args.surface_at2m)

