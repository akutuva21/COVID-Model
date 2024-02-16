import re

def if_else_to_bnet(file_path):
    # Step 1: Define a function that takes in the boolean network as a string.
    conditions_regex = r"(if|elif)\s+(.+)\s*:"
    assignments_regex = r"^\s+(.+)\s*=\s*(.+)$"
    
    # Step 2: Use regular expressions to extract the conditions and assignments from the if/else/elif statements.
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    conditions = []
    assignments = []
    for line in lines:
        if re.match(conditions_regex, line):
            conditions.append(re.match(conditions_regex, line).group(2))
        elif re.match(assignments_regex, line):
            assignments.append(re.match(assignments_regex, line).groups())
    
    # Step 3: Create a dictionary to store the variables and their boolean functions.
    variables = {}
    for assignment in assignments:
        variables[assignment[0]] = assignment[1]
    
    # Step 4: Iterate through the extracted conditions and assignments and update the dictionary accordingly.
    for i in range(len(conditions)):
        condition = conditions[i]
        condition = re.sub(r"\b(and|or|not)\b", r" \1 ", condition)
        condition = re.sub(r"\s+", " ", condition).strip()
        print(f"Condition {i}: {condition}")
        print(f"Variables: {variables}")
        # Evaluate the condition
        if i == 0:
            # If the first condition is True, assign the variable.
            if eval(condition, variables):
                variable = re.search(r"\b([a-zA-Z0-9_]+)\s*=", assignments[i][0]).group(1)
                variables[variable] = eval(assignments[i][1], variables)
        else:
            # Otherwise, iterate through the previous conditions and assign the variable with the first True condition.
            found_true_condition = False
            for j in range(i):
                if eval(conditions[j], variables):
                    variable = re.search(r"\b([a-zA-Z0-9_]+)\s*=", assignments[j][0]).group(1)
                    variables[variable] = eval(assignments[j][1], variables)
                    found_true_condition = True
                    break
            # If none of the previous conditions were True, assign the variable with the current condition.
            if not found_true_condition:
                variable = re.search(r"\b([a-zA-Z0-9_]+)\s*=", assignments[i][0]).group(1)
                variables[variable] = eval(assignments[i][1], variables)
    
    # Step 5: Convert the dictionary to .bnet form and return the result.
    bnet_output = ""
    for variable, function in variables.items():
        bnet_output += variable + ", " + str(function).replace(" and ", " & ").replace(" or ", " | ").replace("not ", "!") + "\n"
    return bnet_output.strip()

result = if_else_to_bnet("code_file.txt")
print(result)