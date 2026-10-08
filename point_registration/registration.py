# LMRG Study 4 Registration

import pathlib
import re
import time
from collections import namedtuple
from pathlib import Path
from typing import Hashable

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial import distance_matrix
from sklearn.neighbors import KDTree

# ## Helper functions
# Define helper functions for extracting scale, angle, and translation from the registration results.

def get_scale(x) -> pd.Series:
    if x[1] is None:
        return pd.Series(np.nan)
    else:
        return pd.Series(x[1][0])


def get_angle(x) -> pd.Series:
    if x[1] is None:
        return pd.Series(np.nan)
    else:
        r = x[1][1]
        return pd.Series(
                np.arctan2(r[1, 0], r[0, 0])
        )


def get_translation(x) -> pd.Series:
    if x[1] is None:
        return pd.Series(np.array([np.nan, np.nan], dtype=np.float64))
    else:
        return pd.Series(x[1][2])


# ## Map file names
# Map the ground truth file names to the friendly file names used in the released images.

friendly_names = {
    'celegans_dyn-90_ceff-0_label.ics.ome.tiff': 'fish1',
    'celegans_dyn-90_ceff-90_label.ics.ome.tiff': 'fish2',
    'celegans_dyn-10_ceff-0_label.ics.ome.tiff': 'fish3',
    'celegans_dyn-10_ceff-90_label.ics.ome.tiff': 'fish4',

    'Q10.13': 'fish1',
    'Q11.13': 'fish2',
    'Q12.13': 'fish3',
    'Q13.13': 'fish4',

    'out_c00_dr10_label.tif': 'nuclei3',
    'out_c00_dr90_label.tif': 'nuclei1',
    'out_c90_dr10_label.tif': 'nuclei4',
    'out_c90_dr90_label.tif': 'nuclei2',

    'Q6.13': 'nuclei1',
    'Q7.13': 'nuclei2',
    'Q8.13': 'nuclei3',
    'Q9.13': 'nuclei4',

    }

PSF_MAP_XY = {  # micrometers, max measured (not theoretical)
    'fish': 0.290,
    'nuclei': None,
}
PSF_MAP_Z = {  # micrometers, measured (not theoretical)
    'fish': 0.182,
    'nuclei': None,
}

# ## Load data
# Load the ground truth data.

ground_truth_coords = pd.read_csv(
    'ground_truth/ground_truth_coords_scale_corrected.csv',
    index_col=0,
)
ground_truth_coords['ground_truth_name'] = (
    ground_truth_coords['path']
    .apply(lambda filename: friendly_names[pathlib.Path(filename).name])
)

# Load the test (submission) data and assigned it the ground truth by its friendly name.
raw_data = pd.read_csv('./all_data_deidentified.csv').dropna(subset=['x', 'y', 'z'])
corrected_data = pd.read_csv('./all_data_deidentified_scale_corrected.csv')
registered_data = pd.read_csv('./all_data_deidentified_scale_corrected_with_registration.csv')
for test_data in [raw_data, corrected_data, registered_data]:
    test_data['filename'] = test_data['csv_path'].apply(lambda f: Path(f).name)
    for f, d in test_data.groupby('csv_path'):
        regex = re.match(r'.*((?:fish|nuclei)[1-4]).*', str(f).lower())
        if regex is not None:
            name = regex.group(1)
            test_data.loc[
                test_data.loc[:, 'csv_path'] == f, 'ground_truth_name'
            ] = name
        elif 'fish' in str(f).lower():
            regex = re.match(r'.*(Q1[0-4]\.13).*', str(f))
            if regex is not None:
                name = friendly_names[regex.group(1)]
                test_data.loc[
                    test_data.loc[:, 'csv_path'] == f, 'ground_truth_name'
                ] = name
        elif 'nuclei' in str(f).lower():
            regex = re.match(r'.*(Q[6-9]\.13).*', str(f))
            if regex is not None:
                name = friendly_names[regex.group(1)]
                test_data.loc[
                    test_data.loc[:, 'csv_path'] == f, 'ground_truth_name'
                ] = name
        else:
            print('No match:', f)

    # Ensure that each `ground_truth_name` has been assigned, that is, no `ground_truth_name` is empty. `True` if correct
    for f in test_data['csv_path'].unique():
        if 'fish' not in f.lower():
            continue
        try:
            td_orig = test_data[test_data.loc[:, 'csv_path'] == f]
            gt_orig = ground_truth_coords[
                ground_truth_coords['ground_truth_name'] == td_orig['ground_truth_name'].iloc[0]
            ][['x', 'y', 'z']].values.astype('float32')
            td_orig = td_orig[['x', 'y', 'z']].values.astype('float32')
            if (re_match := re.match(r'(R_[0-9A-Za-z]+)_.*', pathlib.Path(f).name)) is not None:
                print(re_match.group(1))
        except ValueError as ve:
            print('ValueError:', ve)

# ## Registration (only in x-y)
# Perform rigid registration (with scale) on all the data sets, then extract the parameters.

results = {'result': [], 'filename': []}
scales_map_xy = {  # micrometers
    'fish': 0.1616,
    'nuclei': 0.124,
}
scales_map_z = {  # micrometers
    'fish': 0.200,
    'nuclei': 0.200,
}

RegistrationResults = namedtuple('RegistrationResults', [
    'test_coords_tf',
    'test_tf',
    'test_coords',
    'ground_truth_coords',
    'outliers',
])


def register(args: tuple[Hashable, pd.DataFrame]) -> 'dict[str, str | RegistrationResults]':
    f, _ = args

    results = {}
    results['filename'] = pathlib.Path(str(f)).name

    td_orig = test_data[test_data.loc[:, 'csv_path'] == f]
    gt_orig = ground_truth_coords[
        ground_truth_coords['ground_truth_name'] == td_orig['ground_truth_name'].iloc[0]
    ][['x', 'y', 'z']].values.astype('float64')
    gt = gt_orig.copy()

    results['result'] = RegistrationResults(**{
        'test_coords_tf': td_orig.loc[:, ['x', 'y', 'z']].dropna().copy(deep=True).values,
        'test_tf': None,
        'test_coords': td_orig.loc[:, ['x', 'y', 'z']].dropna().values,
        'ground_truth_coords': gt,
        'outliers': {'0': []},
    })
    return results
    

t0 = time.time()
results_map = map(
    register,
    test_data.groupby('csv_path'),
)
results_df = pd.DataFrame(columns=['result', 'filename'])
for i, r in enumerate(results_map):
    results_df.loc[i, 'result'] = r['result']
    results_df.loc[i, 'filename'] = r['filename']

# Perform linear sum assignment on the registration results, with test points sorted by increasing distance to nearest neigbor in the ground truth. In the case where the scale is close enough to 1.0, don't use the transformed test coordinates and instead use the original submitted coordinates.
t0 = time.time()
distances = {}
lsa_indices = {}
# for filename in results_df.dropna(subset='scale')['filename']:
for filename in results_df['filename']:
    test_coords_tf, test_tf, test_coords, gt, outliers = results_df.query(
        'filename == @filename'
    )['result'].iloc[0]
    kd = KDTree(gt)
    td = test_coords
    pattern = re.match(r'.*(fish|nuclei).*', filename.lower())
    if pattern is not None:
        scale_xy = scales_map_xy[pattern.group(1)]
    else:
        continue
    distance, index = kd.query(test_coords_tf, k=1, return_distance=True)
    distances[filename] = distance
    dm = distance_matrix(gt, test_coords_tf)
    lsa_indices[filename] = {}
    lsa_indices[filename]['dm'] = dm
    lsa_indices[filename]['lsa'] = linear_sum_assignment(dm)

results_df.iloc[[0], :]['result'].item()

results_df_success = results_df.drop(results_df[results_df['result'].apply(lambda x: x.test_coords_tf is None)].index)
results_df_success = pd.concat(
    [
        results_df_success, 
        results_df_success['filename'].str.extract(r".*(?P<dataset>(?P<category>fish|nuclei)[1-4]).*"),
    ],
    axis=1,
)

ground_truth = {
    category:
    {
        dataset: data for (dataset, _), data in ground_truth_coords.groupby(['ground_truth_name', 'path']) if category in dataset  # ty: ignore
    }
    for category in ['fish', 'nuclei']
}


def check_if_close(x):
    # Original ground truth
    ratio = (ground_truth[x['category'].item()][x['dataset'].item()].values[:, 2].std() / 
     # Original z-coordinates (all)
     np.vstack(
         [
             x['result'].item().test_coords,
             *list(x['result'].item().outliers.values()),
         ],
     )[:, 2].std())
    # z scale in micrometers
    return pd.Series(
        [
            ratio,
            np.isclose(ratio, scales_map_z[x['category'].item()], rtol=0.1, atol=0.02),
            np.isclose(ratio, 1.0, rtol=0.1, atol=0.02),
            np.isclose(ratio, scales_map_z[x['category'].item()] / scales_map_xy[x['category'].item()], rtol=0.1, atol=0.02),
        ]
    )


def do_lsa(row):
    _, x = row
    series = pd.DataFrame.from_dict({'lsa': [linear_sum_assignment(
        distance_matrix(
            x['result'].ground_truth_coords, 
            x['result'].test_coords
        )
    )]})
    df = pd.concat([pd.DataFrame(np.array([x.values])), series], axis=1, ignore_index=True)
    df.columns = [*x.index, *series.columns]
    return df


lsa = pd.concat(
    map(do_lsa, results_df_success.dropna().query('category == "fish"').iterrows()),
    ignore_index=True,
)

# Raw data. This data has **not** been registered.
def gen_raw_data(data) -> pd.DataFrame:
    raw_data_list = []
    for filename in data['filename'].unique():
        df = pd.DataFrame(columns=['filename', 'result', 'ground_truth_name'])
        test = data[data['filename'] == filename]
        df['filename'] = [filename]
        df['ground_truth_name'] = [test['ground_truth_name'].iloc[0]]
        ground = ground_truth_coords[ground_truth_coords['ground_truth_name'] == df['ground_truth_name'].iloc[0]]
        df['result'] = [RegistrationResults(
            test_coords_tf=test[['x', 'y', 'z']].dropna(axis=0, how='any').values,
            test_tf=None,
            test_coords=test[['x', 'y', 'z']].dropna(axis=0, how='any').values,
            ground_truth_coords=ground[['x', 'y', 'z']].dropna(axis=0, how='any').values,
            outliers={'0': []},
        )]
        raw_data_list.append(df)
    return pd.concat(raw_data_list, ignore_index=True)
raw_data = gen_raw_data(raw_data)
raw_data['analysis_level'] = 'raw_data'


def iterrows_preserve_dtypes(df):
    return ((i, df.loc[[i], :]) for i in df.index)


def nn_dist(df):
    """Returns mean and std of nn distribution."""
    means = []
    stdevs = []
    for _, row in iterrows_preserve_dtypes(df):
        try:
            ground = row['result'].item().ground_truth_coords
            test = row['result'].item().test_coords
            kd = KDTree(ground)
            distance, _ = kd.query(test, k=1, return_distance=True)
            means.append(np.mean(distance))
            stdevs.append(np.std(distance))
        except ValueError as ve:
            print('ValueError (l507)', ve, row['filename'].item())
            means.append(np.nan)
            stdevs.append(np.nan)
    df_out = df.copy(deep=True)
    df_out['nn_mean'] = means
    df_out['nn_std'] = stdevs
    return df_out.drop(columns=['result', 'ground_truth_name', 'analysis_level']).copy(deep=True)


def lsa_dist_and_jaccard(df):
    """Returns mean and std of lsa distribution."""
    dfs = map(_lsa_dist_and_jaccard_helper, iterrows_preserve_dtypes(df))

    df_out = pd.concat(dfs)
    return df_out


def _lsa_dist_and_jaccard_helper(df_iter):
    """Helper function to enable data parallelism in `lsa_dist_and_jaccard`."""
    _, row = df_iter
    row = row.copy(deep=True)
    try:
        ground = row['result'].item().ground_truth_coords.astype('float64')
        test = row['result'].item().test_coords_tf.astype('float64')
        dm = distance_matrix(ground, test)
        lsa = linear_sum_assignment(dm)
        displacement = ground[lsa[0]] - test[lsa[1]]
        distance = np.sqrt(np.sum(displacement**2, axis=1))

        lsa_mean = np.mean(distance)
        lsa_std = np.std(distance)
        lsa_mse = np.mean(distance**2)
        row['lsa'] = [lsa]
        row['lsa_mean'] = lsa_mean
        row['lsa_std'] = lsa_std
        row['lsa_mse'] = lsa_mse

        category_match = re.match(r'(fish|nuclei)[1-4]', row['ground_truth_name'].item())
        if category_match is not None:
            category = category_match.group(1)
            j_distance = np.sqrt(np.sum((displacement / [[PSF_MAP_XY[category]] * 2 + [PSF_MAP_Z[category]]])**2, axis=1))
            tp = np.sum(j_distance <= 1)
            fp = (test.shape[0] - ground.shape[0] if test.shape[0] >= ground.shape[0] else 0) \
               + sum(map(len, row['result'].item().outliers.values())) \
               + np.sum(j_distance > 1)
            fn = ground.shape[0] - test.shape[0] if test.shape[0] < ground.shape[0] else 0
            jac = tp / (tp + fp + fn)
            if jac < 0.5:
                print(row['filename'].item())
            row['tp'] = [tp]
            row['fp'] = [fp]
            row['fn'] = [fn]
            row['jac'] = [jac]
        else:
            raise TypeError

    except (AttributeError, TypeError, ValueError):
        if 'lsa' not in locals():
            row['lsa'] = [None]

        for var in ['lsa_mean', 'lsa_std', 'lsa_mse', 'tp', 'fp', 'fn', 'jac']:
            if var in locals():
                row[var] = locals()[var]
            else:
                row[var] = [np.nan]

    return row.drop(columns=['result', 'ground_truth_name', 'analysis_level']).copy(deep=True)


jac = lsa_dist_and_jaccard(raw_data)

nn = nn_dist(raw_data)

raw_data_results = raw_data.join(jac.set_index('filename'), on='filename').join(nn.set_index('filename'), on='filename')
raw_data_results['scale_xy'] = raw_data_results['result'].apply(get_scale)
raw_data_results[['translation_x', 'translation_y']] = raw_data_results['result'].apply(get_translation)
raw_data_results['angle_xy'] = raw_data_results['result'].apply(get_angle)

# Add um units
for units, columns in zip(
        ['um', 'um^2', 'px'],
        [
            ['lsa_mean', 'lsa_std', 'nn_mean', 'nn_std'],
            ['lsa_mse'],
        ],
):
    for col in columns:
        if col in raw_data_results.columns:
            raw_data_results[col + ' (' + units + ')'] = raw_data_results[col].copy()
            raw_data_results.drop(col, axis=1, inplace=True)

corrected_data = gen_raw_data(corrected_data)
corrected_data['analysis_level'] = 'corrected'

fully_transformed_data = gen_raw_data(registered_data)
fully_transformed_data['analysis_level'] = 'fully_transformed'

nn_corrected = nn_dist(corrected_data)
lsa_corrected = lsa_dist_and_jaccard(corrected_data)

corrected_results = corrected_data.join(lsa_corrected.set_index('filename'), on='filename').join(nn_corrected.set_index('filename'), on='filename')
corrected_results['scale_xy'] = corrected_results['result'].apply(get_scale)
corrected_results[['translation_x', 'translation_y']] = corrected_results['result'].apply(get_translation)
corrected_results['angle_xy'] = corrected_results['result'].apply(get_angle)

nn_full = nn_dist(fully_transformed_data)
lsa_full = lsa_dist_and_jaccard(fully_transformed_data)

fully_transformed_results = fully_transformed_data.join(lsa_full.set_index('filename'), on='filename').join(nn_full.set_index('filename'), on='filename')
fully_transformed_results['scale_xy'] = fully_transformed_results['result'].apply(get_scale)
fully_transformed_results[['translation_x', 'translation_y']] = fully_transformed_results['result'].apply(get_translation)
fully_transformed_results['angle_xy'] = fully_transformed_results['result'].apply(get_angle)

# Add units to column values that should have real units
for units, columns in zip(
        ['um', 'um^2', 'px'],
        [
            ['lsa_mean', 'lsa_std', 'nn_mean', 'nn_std'],
            ['lsa_mse'],
        ],
):
    for col in columns:
        if col in corrected_results.columns:
            corrected_results[f'{col} ({units})'] = corrected_results[col].copy()
            corrected_results.drop(col, axis=1, inplace=True)

        if col in fully_transformed_results.columns:
            fully_transformed_results[f'{col} ({units})'] = fully_transformed_results[col].copy()
            fully_transformed_results.drop(col, axis=1, inplace=True)

print(corrected_results[raw_data_results['lsa_mse (um^2)'] == corrected_results['lsa_mse (um^2)']].query('filename.str.contains("fish|FISH")'))

output_stats = pd.concat(
    [raw_data_results, corrected_results, fully_transformed_results],
    ignore_index=True,
)
output_stats['id'] = output_stats['filename'].str.extract(r'(R_[0-9A-Za-z]+).*')

# EXPORTED_RESULTS
output_stats.drop(columns=['filename', 'result', 'lsa']).to_csv(
    'registration_stats.csv', index=False)
