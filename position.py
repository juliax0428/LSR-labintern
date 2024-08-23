import pandas as pd

def convert_column_type(column):
    """
    Attempt to convert a column to float.
    If it fails, convert to string.
    """
    try:
        # Try to convert the column to float
        return column.astype(float)
    except ValueError:
        try:
            # If it fails, convert the column to string
            return column.astype(str)
        except Exception as e:
            print(f"Error converting column: {e}")
            return column

def process_csv(file_path):
    """
    Load the CSV file, fill all NaN values with 0, and attempt to convert each column.
    Then, save the updated DataFrame back to a CSV file.
    """
    # Load the CSV file into a DataFrame
    df = pd.read_csv(file_path, delimiter=';',na_values=" ").fillna(value = 0)
    
    # Convert each column's data type
    for column_name in df.columns:
        df[column_name] = convert_column_type(df[column_name])
    
    # Save the updated DataFrame back to a CSV file
    output_file = file_path.replace('.csv', '_converted.csv')
    df.to_csv(output_file, index=False, sep=';')
    print(f"Processed file saved as: {output_file}")

# Example usage
file_path = '/Users/xxz/Desktop/LSR-labintern/czernik20.csv'  # Update this with the correct path
process_csv(file_path)