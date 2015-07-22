### Script to extract all gene features

# We split the processing so we can do it un parallel

# extract from 0 to 1000
#run extract_features-run0.py

# extract from 1000 to 2000
#run extract_features-run1000.py

# extract from 2000 to 3000
#run extract_features-run2000.py

# extract from 3000 to the end
#run extract_features-run3000.py

# We load the 4 files
import pandas as pd

df_features0      = pd.read_csv('./processed/features_genes0.csv', sep=' ')
df_features1000      = pd.read_csv('./processed/features_genes1000.csv', sep=' ')
df_features2000      = pd.read_csv('./processed/features_genes2000.csv', sep=' ')
df_features3000      = pd.read_csv('./processed/features_genes3000.csv', sep=' ')

df_features = df_features0
df_features = df_features.append(df_features1000)
df_features = df_features.append(df_features2000)
df_features = df_features.append(df_features3000)
