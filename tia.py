def main():
    print("Hello from glam!")


if __name__ == "__main__":
    main()

# Importing

from pyteomics import parser, mass
import pandas as pd
import csv



# Load digested peptides from a CSV file.
# df is the data frame


def load_peptides(peptide_file):
    try:
        peptides_df = pd.read_csv(peptide_file)
        if 'peptide' not in peptides_df.columns or 'monoisotopic mass' not in peptides_df.columns:
            raise ValueError("The peptide file must contain 'peptide' and 'monoisotopic mass' columns.")
        return peptides_df.to_dict('records')
    except Exception as e:
        print(f"Error loading peptide file: {e}")
        exit(1)

# Load glycans from the user inputted glycan CSV file
def load_glycans(glycan_file):
    try:
        glycans_df = pd.read_csv(glycan_file)
        if 'glycan' not in glycans_df.columns or 'monoisotopic mass' not in glycans_df.columns:
            raise ValueError("The glycan file must contain 'glycan' and 'monoisoptic mass' columns.")
        return glycans_df.to_dict('records')
    except Exception as e:
        print(f"Error loading glycan file: {e}")
        exit(1)

# Might not need this depending on what brooks has done but this function is to calculate peptide mass using Pyteomics. 
# I wrote this without turning it into a variable so can change for more claarity if needed

def calculate_peptide_mass(peptide):
    try:
        return mass.fast_mass(peptide)
    except Exception as e:
        print(f"Error calculating mass for peptide '{peptide}': {e}")
        return None

# Function to generate glycopeptides and calculate their monoisotopic mass
# Taking into consideration when a glycan attaches to a peptide a covalent bond is formed.
# The bond releases H₂O so reducing the overall mass by 18.01528 Da
# Should I write this as a try/except function? i.e for errors


def generate_glycopeptides(peptides, glycans):
    glycopeptides = []
    water_mass = 18.01528
    for peptide in peptides:
        peptide_mass = calculate_peptide_mass(peptide)
        if peptide_mass is None:
            continue
        for glycan in glycans:
            glycopeptide_sequence = peptide + "+" + glycan['glycan']
            glycopeptide_mass = peptide_mass + glycan['monoisotopic mass'] - water_mass
            glycopeptides.append({
                'glycopeptide': glycopeptide_sequence,
                'monoisotopic mass': glycopeptide_mass
            })
    return glycopeptides

# Main function
def main():
    # Load peptides 
    peptide_file = input("Enter the path to the peptide file: ").strip()
    peptides = load_peptides(peptide_file)

    # Load glycans 
    glycan_file = input("Enter the path to the glycan file: ").strip()
    glycans = load_glycans(glycan_file)

    # Generate glycopeptides
    glycopeptides = generate_glycopeptides(peptides, glycans)

    # Save the output to a CSV file
    output_df = pd.DataFrame(glycopeptides)
    output_file = "glycopeptides_output.csv"
    output_df.to_csv(output_file, columns=['glycopeptide', 'monoisotopic mass'], index=False)
    print(f"Glycopeptides saved to {output_file}")

if __name__ == "__main__":
    main()


