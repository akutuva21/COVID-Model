import argparse


def parse_if_else_statements(code_str):
    if_else_statements = {}
    for line in code_str.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith(('if', 'elif')):
            var_name = line.split('==')[1].strip().replace('"', '').replace("'", "")
            if_else_statements[var_name] = []
        elif line.startswith('else'):
            if_else_statements['else'] = []
        elif line.startswith('x['):
            if_else_statements[var_name].append(line.strip())
        else:
            if_else_statements[var_name][-1] += line.strip()
    return if_else_statements


def generate_bnet_file(if_else_statements, output_file):
    with open(output_file, 'w') as f:
        f.write('Variables:\n')
        for var_name, expr_list in if_else_statements.items():
            f.write(f"{var_name}, {' & '.join(expr_list)}\n")
        f.write('\nInteractions:\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate boolean network files from if-else statements.')
    parser.add_argument('input_file', type=str, help='path to the input file')
    parser.add_argument('output_file', type=str, help='path to the output file')
    args = parser.parse_args()

    with open(args.input_file, 'r') as f:
        code_str = f.read()

    if_else_statements = parse_if_else_statements(code_str)
    generate_bnet_file(if_else_statements, args.output_file)