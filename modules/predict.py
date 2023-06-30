
import numpy as np
import pandas as pd
import numpy as np
import copy
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.compose import make_column_selector as selector
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer



def get_spatial_fields(gdf):
    
    spatial_gdf = pd.DataFrame()
    centroid = gdf['geometry'].centroid
    spatial_gdf['area'] = gdf['geometry'].area
    spatial_gdf['perimeter'] = gdf['geometry'].length
    spatial_gdf['x'] = centroid.x
    spatial_gdf['y'] = centroid.y
    return spatial_gdf


def create_train_test(complete_data, target_name, test_size=0.4):
    """
    The function splits the data to in training and test set according to a specific test size.
    """
    y1 = complete_data[target_name]
    X1 = complete_data.drop(columns=[target_name])

    X_train, X_test, y_train, y_test = train_test_split(
            X1, y1, test_size, random_state=42
        )
    return X_train, X_test, y_train, y_test

def fill_missing_values(target_name, gdf, selected_columns, clf=RandomForestClassifier(max_depth=5, n_estimators=200)):
    """
    The function takes a geodataframe, the name of the target feature, the list of selected columns and the classifier 
    to use to fill the missing values of the selected column.
    Returns the selected columns with the missing values filled.
    """
    data = copy.deepcopy(gdf)
    spatial_fields = get_spatial_fields(data)
    for field in spatial_fields:
        data[field] = spatial_fields[field]
    
    if target_name in data.columns:     # select the columns that are in selected columns and the target feature
        data = data[selected_columns]
        # create a column containing the index

    else:
        data = data[selected_columns + [target_name]]

    # keep track of the original order of the columns

    complete_data = copy.deepcopy(data[data[target_name].notna()])

    incomplete_data = copy.deepcopy(data[data[target_name].isna()])
    y = complete_data[target_name]
    X = complete_data.drop(columns=[target_name])

    X_pred = incomplete_data.drop(columns=[target_name])

    imp = SimpleImputer(missing_values=np.nan, strategy='mean')

    numerical_columns_selector = selector(dtype_exclude=object)
    categorical_columns_selector = selector(dtype_include=object)

    numerical_columns = numerical_columns_selector(X)
    categorical_columns = categorical_columns_selector(X)

    for numerical_col in numerical_columns: # for numerical columns missing values have to be replaced (while for categorical columns we can keep nan values!)
        if X_pred[numerical_col].isna().any():
            imp.fit(X_pred[X_pred[numerical_col].notna()][numerical_col].values.reshape(-1,1)) 
            X_pred.loc[X_pred[numerical_col].isna(), numerical_col] = imp.transform(X_pred[X_pred[numerical_col].isna()][numerical_col].values.reshape(-1,1)).reshape(-1,)
        # First, let’s create the preprocessors for the numerical and categorical parts.
    
    categorical_preprocessor = OneHotEncoder(handle_unknown="ignore")
    numerical_preprocessor = StandardScaler()

    # Now, we create the transformer and associate each of these preprocessors with their respective columns.
    preprocessor = ColumnTransformer([
        ('one-hot-encoder', categorical_preprocessor, categorical_columns),
        ('standard_scaler', numerical_preprocessor, numerical_columns)],
        remainder="passthrough")
    clf = make_pipeline(preprocessor, clf)
    clf.fit(X, y)
    score = clf.score(X, y)
    model = str(clf.named_steps['randomforestclassifier'])
    model_and_score = str("Model: %s" % model) + str(" Score: %.4f" % score)

    #print(X_pred.isna().values)
    y_pred= clf.predict(X_pred)
    # put y_pred and X_pred in the same dataframe 
    incomplete_data[target_name] = y_pred
    incomplete_data[target_name + '_IsPredicted'] = 'predicted'
    incomplete_data[target_name + '_comments'] = model_and_score
    complete_data[target_name + '_IsPredicted'] = 'user'
    #add the incomplete data to the complete data
    complete_data = pd.concat([complete_data, incomplete_data], axis=0)
    complete_data.sort_index(inplace=True)
    completed_column = complete_data[target_name]
    is_predicted_column = complete_data[target_name + '_IsPredicted']
    comments_column = complete_data[target_name + '_comments']
    return completed_column, is_predicted_column, comments_column

