import pathmodel

# pathmodel infer -i data.lp -o output_folder -s 100

INPUT_FILE = ''
OUTPUT_DIR = '../Files/Outputs/oxylipins'

pathmodel.pathmodel_analysis(INPUT_FILE, OUTPUT_DIR, step_limit=100)
# pathmodel_plot -i OUT