'''
NAME
    metabolite_production

VERSION
    1.0
    
TYPE
    Homework

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/metabolite_production.py

DESCRIPTION
    Returns the genes with the highest values of each category in a dataframe
    about a metabolite production in specific conditions.

CATEGORY
    Pandas
    
REQUIREMENTS
    - Software - 
        You must have Python and BioPython installed
    - Repository structure -
        For using this script successfully, you should have the 
        following repository structure in front of you:
        |- data      # It contains data to be used as input
        |- docs      # It contains documentation
        |- results   # Output files will be created in this directory
        |- src       # The script must be executed from this directory

USAGE
    python src/metabolite_production.py

ARGUMENTS
    None
    
SEE ALSO
    usepand.py
'''

# Import libraries
import pandas as pd

# Create costs column for the dataframe
costs = pd.Series([ 3.5, 5, 7, 4.3],
                      index=['gene1', 'gene2', 'gene3','gene5'],
                      name='costs')

# Create production at 30 degrees column for the dataframe
production_30 = pd.Series([5, 11, 4, 7, 2],
                        index=["gene1", "gene2", "gene3", "gene4", "gene5"],
                        name="production at 35°C")

# Create production at 35 degrees column for the dataframe
production_35 = pd.Series([3, 7, 9, 4, 6],
                        index=["gene1", "gene2", "gene3", "gene4", "gene5"],
                        name="production at 35°C")

# Create cost-benefit production
cost_benefit = pd.DataFrame({'costs': costs,
                       'production_at_30': production_30, 
                       'production_at_35': production_35})

# Calculate unitary costs
unitary_cost_30 = cost_benefit.costs / cost_benefit.production_at_30
unitary_cost_35 = cost_benefit.costs / cost_benefit.production_at_35

# Add unitary costs columns to the dataframe
cost_benefit.insert(2, "unitary_cost_at_30", unitary_cost_30)
cost_benefit.insert(4, "unitary_cost_at_35", unitary_cost_35)

# Print results
print(cost_benefit.idxmax())
