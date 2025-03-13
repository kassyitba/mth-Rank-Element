import random
import statistics


class Custodian:
    def __init__(self, data, id):
        """
        Initialize the custodian with its dataset.
        - data: A list of numbers in their original order.
        - id: An identifier for the custodian.

        The constructor stores:
          - raw_data: The original dataset.
          - data: The sorted dataset (R_j) for internal operations.
          - unique_data: The sorted set of unique elements (U_j) for local median calculations.
        """
        self.id = id
        self.raw_data = data[:]  # Original dataset for display.
        self.data = sorted(data)  # Sorted dataset R_j for internal use.
        self.unique_data = sorted(set(data))  # Unique sorted dataset U_j.

    def display_original(self):
        """Display the original dataset for this custodian."""
        print(f"Custodian {self.id} - Original Dataset: {self.raw_data}")

    def display_sorted(self):
        """Display the sorted dataset for this custodian."""
        print(f"Custodian {self.id} - Sorted Dataset: {self.data}")

    def get_local_median_in_range(self, x_lower, x_upper):
        """
        Compute the local median m_j^i over unique elements in the range (x_lower, x_upper].
        If there are such elements, return the median and flag (c_j^i = 1).
        Otherwise, return 0 and flag 0.
        """
        # Select unique elements within the range (x_lower, x_upper]
        candidates = [val for val in self.unique_data if x_lower < val <= x_upper]
        if candidates:
            median_val = statistics.median(candidates)
            return median_val, 1  # Valid contribution.
        else:
            return 0, 0  # No valid candidate in range.

    def count_elements_leq(self, pivot):
        """
        Count the number of elements in the sorted dataset (R_j) that are <= pivot.
        Since R_j is sorted, the loop stops as soon as an element exceeds the pivot.
        """
        count = 0
        for val in self.data:
            if val <= pivot:
                count += 1
            else:
                break
        return count


def simulate_protocol(custodians, target_rank, max_iterations=1000, tolerance=1e-6):
    """
    Simulate the distributed selection protocol to find the m-th smallest element.

    Parameters:
      - custodians: List of Custodian objects.
      - target_rank: The target m-th smallest element (1-indexed).
      - max_iterations: Maximum iterations allowed.
      - tolerance: Convergence threshold for the search interval.

    The protocol iteratively:
      1. Establishes a global search range from the minimum and maximum values.
      2. Each custodian computes a local median for values within the current range.
      3. A global pivot is computed via simulated secure summation.
      4. Each custodian counts the number of elements <= the pivot.
      5. The search range is updated based on the aggregated count relative to target_rank.

    Returns the approximate m-th smallest element.
    """
    # Establish the initial global search bounds from all custodians.
    global_min = min(custodian.data[0] for custodian in custodians if custodian.data)
    global_max = max(custodian.data[-1] for custodian in custodians if custodian.data)
    x_lower = global_min
    x_upper = global_max
    iteration = 0

    while iteration < max_iterations:
        iteration += 1
        print(f"\nIteration {iteration}:")
        print(f"  Current search range: ({x_lower}, {x_upper}]")

        # Each custodian computes its local median for values within the current range.
        local_medians = []
        contribution_flags = []
        for custodian in custodians:
            m_j, c_j = custodian.get_local_median_in_range(x_lower, x_upper)
            local_medians.append(m_j)
            contribution_flags.append(c_j)
            #print(f"  Custodian {custodian.id}: Local median in range = {m_j} with flag = {c_j}")
            print(f"  Custodian {custodian.id}: Local median in range = {m_j} ")

        # Global pivot calculation via simulated secure summation.
        total_median_sum = sum(local_medians)
        total_contribution = sum(contribution_flags)
        if total_contribution == 0:
            # If no valid medians are found, use the average of current bounds.
            P = (x_lower + x_upper) / 2
            print("  No valid local medians found; using fallback pivot (average of bounds).")
        else:
            P = total_median_sum / total_contribution
        print(f"  Global pivot P = {P}")

        # Each custodian counts the number of elements <= P.
        local_counts = []
        for custodian in custodians:
            count = custodian.count_elements_leq(P)
            local_counts.append(count)
            print(f"  Custodian {custodian.id}: Count of elements <= P = {count}")
        global_count = sum(local_counts)
        print(f"  Total count of elements <= P: {global_count}")

        # Check if the count matches the target rank.
        if global_count == target_rank:
            print(f"\nTarget rank reached. The {target_rank}-th smallest element is approximately {P}")
            return P
        elif global_count > target_rank:
            x_upper = P
            print(f"  Global count > target rank. Updating upper bound to P = {P}")
        else:
            x_lower = P
            print(f"  Global count < target rank. Updating lower bound to P = {P}")

        if abs(x_upper - x_lower) < tolerance:
            approx_value = (x_lower + x_upper) / 2
            print("\nConvergence reached based on tolerance.")
            print(f"Approximated {target_rank}-th smallest element is {approx_value}")
            return approx_value

    print("\nMaximum iterations reached without full convergence.")
    approx_value = (x_lower + x_upper) / 2
    print(f"Approximated {target_rank}-th smallest element is {approx_value}")
    return approx_value


# -------------------- Simulation Setup with Random Data -------------------- #

def generate_random_dataset(size, lower_bound=1, upper_bound=100):
    """
    Generate a random dataset of integers.

    Parameters:
      - size: Number of elements in the dataset.
      - lower_bound: Minimum possible integer value.
      - upper_bound: Maximum possible integer value.

    Returns a list of random integers.
    """
    return [random.randint(lower_bound, upper_bound) for _ in range(size)]


# Note: The fixed random seed has been removed so that datasets vary with each run.

# Define the number of custodians.
num_custodians = 3

# Instead of using a fixed size for each custodian, generate varied sizes.
dataset_sizes = [random.randint(5, 15) for _ in range(num_custodians)]

# Create custodians with randomly generated datasets of varied lengths.
custodians = []
for i in range(num_custodians):
    random_data = generate_random_dataset(dataset_sizes[i])
    custodian = Custodian(random_data, id=i + 1)
    custodians.append(custodian)

# Display each custodian's original and sorted datasets.
print("Datasets for Each Custodian:")
for custodian in custodians:
    custodian.display_original()
    custodian.display_sorted()

# Compute and display the combined sorted unique dataset from all custodians.
combined_dataset = []
for custodian in custodians:
    combined_dataset.extend(custodian.raw_data)
combined_unique_sorted = sorted(set(combined_dataset))
print("\nCombined Sorted Unique Dataset from All Custodians:")
print(combined_unique_sorted)

# Define the target rank m (e.g., the 5th smallest element overall).
target_rank = 4

# Run the protocol simulation.
result = simulate_protocol(custodians, target_rank)

print("\nCombined Sorted Unique Dataset from All Custodians:")
print(combined_unique_sorted)

print(f"\nFinal Result: The {target_rank}-th  element is approximately {result}")
