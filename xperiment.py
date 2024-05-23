import subprocess
import csv
import os

# Define the path to the executable and the directory to run it in
executable_path = './test_framework'
working_directory = './build'
#output_filename = 'output.csv'
output_filename = 'output_max.csv'
timeout_duration = 1000  # Timeout in seconds

current_test = 1


# Generate a list for N and k based on growth pattern 10^i + 1
def generate_values(start, end):
    values = []
    i = 0
    value = start
    while value <= end:
        values.append(value)
        i += 1
        value *= 10  # Increasing by a factor of 10 each iteration
    return values

# Run the test framework with varying parameters
def run_tests():
    N_values = generate_values(100000, 100000000)
    k_values = generate_values(100, 100000)
    total_tests = len(N_values) * len(k_values) * 13
    current_test = 1
    with open(output_filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Algorithm', 'N', 'k', 'Alphabet Size', 'Total Time', 'Preprocessing Time'])
        for algo in range(9):  # Algorithm choices from 0 to 12
            for N in N_values:
                for k in k_values:
                    if N == k:  # Check if N is equal to k, and skip this iteration if true
                        print(f"Skipping test where N ({N}) equals k ({k}).")
                        current_test += 1
                        continue
                    # Prepare the subprocess environment
                    env = os.environ.copy()
                    #env['PARLAY_NUM_THREADS'] = '1'
                    if algo in {6, 7, 8}:  # Using a set for faster membership checking
                        env['PARLAY_NUM_THREADS'] = '2'
                    args = [executable_path, str(algo), str(N), str(k), "256"]
                    print(f"Running test {current_test}/{total_tests}: Algorithm {algo}, N {N}, k {k}")
                    try:
                        # Execute the command with the generated arguments and a timeout
                        result = subprocess.run(args, capture_output=True, text=True, check=True, cwd=working_directory, timeout=timeout_duration, env=env)
                        # Extract total time and preprocessing time from the output
                        output = result.stdout.strip().split()
                        total_time, preprocessing_time = output[0], output[1]
                        writer.writerow([algo, N, k, 256, total_time, preprocessing_time])
                    except subprocess.TimeoutExpired:
                        writer.writerow([algo, N, k, 256, 'null', 'null'])
                    except subprocess.CalledProcessError as e:
                        print(f"Error running {args}: {e}")
                        writer.writerow([algo, N, k, 256, 'error', 'error'])
                    current_test += 1

# Main function to start the tests
def main():
    run_tests()
    print(f"Tests completed. Results have been saved to {output_filename}.")

# Execute the script
if __name__ == "__main__":
    main()
