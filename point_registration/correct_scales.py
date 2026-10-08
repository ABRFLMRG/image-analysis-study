from sys import stderr

import numpy as np
import pandas as pd
from pycpd import RigidRegistration


CORRECTIONS: dict[str, dict[str, list[str]]] = {
    'fish': {
        'flip': ['R_3lAJ9xY4kGlL99f'],
        'xyz_pixel_to_um': [
            'R_1q7MrEdSr6yhykH',
            'R_24HIjcCJh6uI3bu',
            'R_3RxOa8kaiyul4bG',
            'R_0638fJAPpmzLtu1',
            'R_31bjqd6Mm8wBxN5',
        ],
        'z_nm_to_um': ['R_2DNRFrAvCDUX1EL'],
        'z_slice_to_um': ['R_2c14tLfUPR1Vnua'],
        'xy_scale_on_z': ['R_1qdHCwPCdNvs9xi'],
        'free_scale_z': ['R_2cCJjlMU7i9XjMQ'],
        'register': [
            'R_1GUZ4XruXifzoPp', 'R_24Nwgngl83ucQ8B', 'R_2DNRFrAvCDUX1EL',
            'R_24HIjcCJh6uI3bu', 'R_3lAJ9xY4kGlL99f', 'R_1q7MrEdSr6yhykH',
            'R_2c14tLfUPR1Vnua', 'R_1M6DoAmYEY3Jrvm', 'R_1ePOgKBqC8pSpzT',
            'R_22s3aTqiX7gbY4u', 'R_2rCMx6wAGE7bFJh', 'R_3j9w1bGwWGd8yOC',
            'R_tY87q7yGKRqww5X', 'R_31bjqd6Mm8wBxN5',
        ],
    },
    'nuclei': {
        'flip': ['R_3lAJ9xY4kGlL99f', 'R_31bjqd6Mm8wBxN5'],
        'xyz_pixel_to_um': ['R_0638fJAPpmzLtu1', 'R_3RxOa8kaiyul4bG'],
        'xyz_4x': ['R_3lAJ9xY4kGlL99f'],
        'free_scale_xy': ['R_2rCMx6wAGE7bFJh', 'R_6E6PxT3N1gr93yh'],
        'free_scale_z': ['R_cHLBH1bftzVtPDH', 'R_6E6PxT3N1gr93yh'],
        'register': ['R_3kpixa7Fm6wlFk2'],
    },
}

SCALES_MAP_XY: dict[str, float] = {  # micrometers
    'fish': 0.1616,
    'nuclei': 0.124,
}
SCALES_MAP_Z: dict[str, float] = {  # micrometers
    'fish': 0.200,
    'nuclei': 0.200,
}

REGISTRATION_PARAMS = {
    'max_iterations': 1_000_000,
}


def flip(input_coords: pd.DataFrame) -> pd.DataFrame:
    coords = input_coords.copy(deep=True)

    coords.loc[:, ['x', 'y']] = pd.concat(
        [coords.loc[:, 'x'].mean() - coords.loc[:, 'x'], coords.loc[:, 'y']], axis=1
    )
    coords.loc[:, ['x', 'y']] = (
        np.array([[0.0, 1.0], [-1.0, 0.0]]) @ coords.loc[:, ['x', 'y']].values.T
    ).T

    return coords


def z_pixel_to_um(input_coords: pd.DataFrame, ground_truth_type: str) -> pd.DataFrame:
    coords = input_coords.copy(deep=True)

    coords.loc[:, ['x', 'y']] *= SCALES_MAP_XY[ground_truth_type]
    coords.loc[:, 'z'] *= SCALES_MAP_Z[ground_truth_type]

    return coords


def xy_pixel_to_um(input_coords: pd.DataFrame, ground_truth_type: str) -> pd.DataFrame:
    coords = input_coords.copy(deep=True)

    coords.loc[:, ['x', 'y']] *= SCALES_MAP_XY[ground_truth_type]

    return coords


def xyz_4x(input_coords: pd.DataFrame) -> pd.DataFrame:
    coords = input_coords.copy(deep=True)

    coords.loc[:, ['x', 'y', 'z']] *= 4.0

    return coords


def z_nm_to_um(input_coords: pd.DataFrame) -> pd.DataFrame:
    coords = input_coords.copy(deep=True)

    if 'R_2DNRFrAvCDUX1EL_fish3' not in input_coords['csv_path'].iloc[0]:
        coords.loc[:, 'z'] /= 1000.0

    return coords


def z_slice_to_um(input_coords: pd.DataFrame, ground_truth_type: str) -> pd.DataFrame:
    coords = input_coords.copy(deep=True)

    coords.loc[:, 'z'] *= SCALES_MAP_Z[ground_truth_type]

    return coords


def xy_scale_on_z(input_coords: pd.DataFrame, ground_truth_type: str) -> pd.DataFrame:
    coords = input_coords.copy(deep=True)

    coords.loc[:, 'z'] *= (
        SCALES_MAP_Z[ground_truth_type] / SCALES_MAP_XY[ground_truth_type]
    )

    return coords


def free_scale_xy(input_coords: pd.DataFrame, ground_truth: pd.DataFrame, ground_truth_type: str) -> pd.DataFrame:
    coords = input_coords.copy(deep=True)

    if 'R_2rCMx6wAGE7bFJh_nuclei1' not in input_coords['csv_path'].iloc[0]:
        scale = (ground_truth[['x', 'y']].std() / input_coords[['x', 'y']].std()).mean()
        coords[['x', 'y']] *= scale

    return coords


def free_scale_z(input_coords: pd.DataFrame, ground_truth: pd.DataFrame, ground_truth_type: str) -> pd.DataFrame:
    coords = input_coords.copy(deep=True)

    if 'R_cHLBH1bftzVtPDH_nuclei4' not in input_coords['csv_path'].iloc[0]:
        scale = ground_truth['z'].std() / input_coords['z'].std()
        coords['z'] *= scale

    return coords


def translate(input_coords: pd.DataFrame, ground_truth: pd.DataFrame, ground_truth_type: str) -> pd.DataFrame:
    coords = input_coords.copy(deep=True)

    translation = ground_truth[['x', 'y', 'z']].mean() - coords[['x', 'y', 'z']].mean()
    coords[['x', 'y', 'z']] += translation

    return coords


def register(input_coords: pd.DataFrame, ground_truth: pd.DataFrame, ground_truth_type: str) -> pd.DataFrame:
    # Copy the input_coords to prevent modifying them.
    coords = input_coords.copy(deep=True)

    # From R_24Nwgngl83ucQ8B only nuclei4 needs registration.
    if 'R_24Nwgngl83ucQ8B_nuclei1' not in input_coords['csv_path'].iloc[0] and \
       'R_24Nwgngl83ucQ8B_nuclei2' not in input_coords['csv_path'].iloc[0] and \
       'R_24Nwgngl83ucQ8B_nuclei3' not in input_coords['csv_path'].iloc[0] and \
       'R_31bjqd6Mm8wBxN5_fish2' not in input_coords['csv_path'].iloc[0]:

        translation = np.array([
            ground_truth[['x', 'y', 'z']].mean(axis=0) - coords[['x', 'y', 'z']].mean(axis=0)
        ])

        try:
            reg = RigidRegistration(
                X=ground_truth.loc[:, ['x', 'y', 'z']].values,
                Y=coords.loc[:, ['x', 'y', 'z']].values,
                s=1.0,  # The data set should have the correct scale at this point.
                t=translation,
                **REGISTRATION_PARAMS,
            )
            coords.loc[:, ['x', 'y', 'z']], _ = reg.register()
        except np.linalg.LinAlgError as err:
            print(input_coords['csv_path'].iloc[0], err)
    else:
        print(f'Skipping {input_coords["csv_path"].iloc[0]}')

    return coords


def main():
    input_data = pd.read_csv('./all_data_deidentified.csv').dropna(subset=['x', 'y', 'z'])
    data = input_data.copy(deep=True)

    ground_truth = pd.read_csv('./ground_truth/ground_truth_coords_scale_corrected.csv')

    data['registered'] = False
    data['needs_translation'] = False

    # Put registered coordinates in this DataFrame.
    registration_coords = data.copy(deep=True)

    for ground_truth_type in ['fish', 'nuclei']:
        for correction_type, ids in CORRECTIONS[ground_truth_type].items():
            for i in ids:
                for n in range(1, 4 + 1):
                    select = data['csv_path'].str.lower().str.contains(
                        f'{i.lower()}_{ground_truth_type}{n}'
                    )
                    # Process only if the data set exists.
                    if select.any():
                        # Corrections are applied one by one. The same
                        # coordinates might receive multiple corrections
                        # but on different iterations.
                        if correction_type == 'flip':
                            data.loc[select, :] = flip(data.loc[select, :])
                            registration_coords.loc[select, :] = data.loc[select, :].copy()
                            data.loc[select, 'needs_translation'] = True
                        elif correction_type == 'xyz_pixel_to_um':
                            data.loc[select, :] = z_pixel_to_um(
                                data.loc[select, :], ground_truth_type
                            )
                            registration_coords.loc[select, :] = data.loc[select, :].copy()
                            data.loc[select, 'needs_translation'] = True
                        elif correction_type == 'xy_pixel_to_um':
                            data.loc[select, :] = xy_pixel_to_um(
                                data.loc[select, :], ground_truth_type
                            )
                            registration_coords.loc[select, :] = data.loc[select, :].copy()
                            data.loc[select, 'needs_translation'] = True
                        elif correction_type == 'xyz_4x':
                            data.loc[select, :] = xyz_4x(data.loc[select, :])
                            registration_coords.loc[select, :] = data.loc[select, :].copy()
                            data.loc[select, 'needs_translation'] = True
                        elif correction_type == 'z_nm_to_um':
                            data.loc[select, :] = z_nm_to_um(data.loc[select, :])
                            registration_coords.loc[select, :] = data.loc[select, :].copy()
                            data.loc[select, 'needs_translation'] = True
                        elif correction_type == 'z_slice_to_um':
                            data.loc[select, :] = z_slice_to_um(
                                data.loc[select, :], ground_truth_type
                            )
                            registration_coords.loc[select, :] = data.loc[select, :].copy()
                            data.loc[select, 'needs_translation'] = True
                        elif correction_type == 'xy_scale_on_z':
                            data.loc[select, :] = xy_scale_on_z(
                                data.loc[select, :], ground_truth_type
                            )
                            registration_coords.loc[select, :] = data.loc[select, :].copy()
                            data.loc[select, 'needs_translation'] = True
                        elif correction_type == 'free_scale_xy':
                            data.loc[select, :] = free_scale_xy(
                                data.loc[select, :],
                                ground_truth.query(f'ground_truth_name == "{ground_truth_type}{n}"'),
                                ground_truth_type,
                            )
                            registration_coords.loc[select, :] = data.loc[select, :].copy()
                            data.loc[select, 'needs_translation'] = True
                        elif correction_type == 'free_scale_z':
                            data.loc[select, :] = free_scale_z(
                                data.loc[select, :],
                                ground_truth.query(f'ground_truth_name == "{ground_truth_type}{n}"'),
                                ground_truth_type,
                            )
                            registration_coords.loc[select, :] = data.loc[select, :].copy()
                            data.loc[select, 'needs_translation'] = True
                        elif correction_type == 'register':
                            registration_coords.loc[select, :] = register(
                                data.loc[select, :].copy(deep=True),
                                ground_truth.query(f'ground_truth_name == "{ground_truth_type}{n}"'),
                                ground_truth_type,
                            )
                            registration_coords.loc[select, 'registered'] = True
                        else:
                            print(
                                f'correction_type {correction_type} not recognized',
                                file=stderr,
                            )
                        print(f'Applied {correction_type} to {i}_{ground_truth_type}{n}')

    # Translate the only the data sets that have been corrected. Since
    # the same data sets can appear more than once in the CORRECTIONS
    # dict, use a 'translated' field to ensure data sets are
    # translated only once, and also to prevent registered data sets
    # from being translated.
    for ground_truth_type in ['fish', 'nuclei']:
        for correction_type, ids in CORRECTIONS[ground_truth_type].items():
            for i in ids:
                for n in range(1, 4 + 1):
                    select = data['csv_path'].str.lower().str.contains(
                        f'{i.lower()}_{ground_truth_type}{n}'
                    )

                    # Translate only if the data exists. Don't
                    # translate R_31bjqd6Mm8wBxN5_fish2 because the
                    # alignment seems better without translation.
                    if select.any() and 'R_31bjqd6Mm8wBxN5_fish2' not in data.loc[select, 'csv_path'].iloc[0]:
                        if data.loc[select, 'needs_translation'].all():
                            data.loc[select, :] = translate(
                                data.loc[select, :],
                                ground_truth.query(f'ground_truth_name == "{ground_truth_type}{n}"'),
                                ground_truth_type,
                            )
                            data.loc[select, 'registered'] = True

    for ground_truth_type in ['fish', 'nuclei']:
        for correction_type, ids in CORRECTIONS[ground_truth_type].items():
            for i in ids:
                for n in range(1, 4 + 1):
                    select = registration_coords['csv_path'].str.lower().str.contains(
                        f'{i.lower()}_{ground_truth_type}{n}'
                    )

                    # Translate only if the data exists. Don't
                    # translate R_31bjqd6Mm8wBxN5_fish2 because the
                    # alignment seems better without translation.
                    if select.any() and 'R_31bjqd6Mm8wBxN5_fish2' not in registration_coords.loc[select, 'csv_path'].iloc[0]:
                        if not registration_coords.loc[select, 'registered'].all():
                            registration_coords.loc[select, :] = translate(
                                registration_coords.loc[select, :],
                                ground_truth.query(f'ground_truth_name == "{ground_truth_type}{n}"'),
                                ground_truth_type,
                            )
                            registration_coords.loc[select, 'registered'] = True

    (
        data
        .drop(['needs_translation', 'registered'], axis=1)
        .to_csv('./all_data_deidentified_scale_corrected.csv', index=False)
    )
    (
        registration_coords
        .drop(['needs_translation', 'registered'], axis=1)
        .to_csv(
            './all_data_deidentified_scale_corrected_with_registration.csv',
            index=False,
        )
    )


if __name__ == '__main__':
    main()
