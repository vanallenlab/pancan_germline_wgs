import sys
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

def main(filename):
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(filename, sep='\t')

    # Define the features (X) and target (y)
    X = df[['risk_hom', 'het', 'safe_hom']]
    y = df['driver_status']

    # Create the logistic regression model
    model = LogisticRegression()

    # Fit the model to the data
    model.fit(X, y)

    # Get the coefficients of the model
    coefficients = model.coef_[0]

    # Calculate the odds ratios (OR)
    odds_ratios = np.exp(coefficients)

    # Create a DataFrame to display results
    results = pd.DataFrame({
        'Feature': X.columns,
        'Coefficient': coefficients,
        'Odds Ratio': odds_ratios
    })

if __name__ == "__main__":
    # Check if the filename argument is provided
    if len(sys.argv) < 2:
        print("Usage: python logistic_regression.py <filename>")
        sys.exit(1)

    # Get the filename from command line argument
    filename = sys.argv[1]

    # Call the main function with the filename argument
    main(filename)
