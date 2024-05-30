import random

def generate_input_file(filename, num_numbers):
    # Open the file for writing
    with open(filename, 'w') as file:
        # Write the number of numbers at the beginning
        file.write(f"{num_numbers} ")
        
        # Generate and write the random numbers
        random_numbers = [random.randint(1, 100) for _ in range(num_numbers)]
        file.write(" ".join(map(str, random_numbers)))
        file.write("\n")

# Specify the filename and number of numbers
filename = "input.txt"
num_numbers = 1000000

# Generate the input file
generate_input_file(filename, num_numbers)
