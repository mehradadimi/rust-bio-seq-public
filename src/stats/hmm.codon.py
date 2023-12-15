from bio import *
import math


class HMM:
    transition_matrix: List[List[float]]
    emission_matrix: List[List[float]]
    initial_state_probs: List[float]
    states: List[int]

    def __init__(self, transition_matrix, emission_matrix, initial_state_probs):
        self.transition_matrix = transition_matrix
        self.emission_matrix = emission_matrix
        self.initial_state_probs = initial_state_probs
        self.states = list(range(len(self.transition_matrix)))

    def viterbi_algorithm(hmm: HMM, obs: List[int]) -> List[int]:
        num_states = len(hmm.states)
        num_obs = len(obs)

        # Initialize Viterbi table and backpointer
        viterbi_table = [[0.0 for _ in range(num_states)] for _ in range(num_obs)]
        backpointer = [[0 for _ in range(num_states)] for _ in range(num_obs)]

        for t in range(num_obs):
            for s in range(num_states):
                if t == 0:
                    viterbi_table[t][s] = (
                        hmm.initial_state_probs[s] * hmm.emission_matrix[s][obs[t]]
                    )

                else:
                    max_prob: float = -1.0
                    max_prev_state = 0

                    for prev_s in range(num_states):
                        # Calculate the transition probability
                        transition_prob = (
                            viterbi_table[t - 1][prev_s]
                            * hmm.transition_matrix[prev_s][s]
                        )

                        if transition_prob > max_prob:
                            max_prob = transition_prob
                            max_prev_state = prev_s

                    # Multiply the max transition probability with the emission probability
                    viterbi_table[t][s] = max_prob * hmm.emission_matrix[s][obs[t]]
                    backpointer[t][s] = max_prev_state

        # Traceback and find best path
        best_path_prob = max(viterbi_table[-1])
        best_path_pointer = viterbi_table[-1].index(best_path_prob)
        best_path = [best_path_pointer]

        for t in range(num_obs - 1, 0, -1):
            best_path.insert(0, backpointer[t][best_path[0]])

        return best_path

    def forward(self, observations: List[int]) -> Tuple[List[List[float]], float]:
        num_states = len(self.states)
        num_obs = len(observations)

        # Logarithmic probabilities
        log_initial_probs = [math.log(p) for p in self.initial_state_probs]
        log_transition_matrix = [
            [math.log(p) for p in row] for row in self.transition_matrix
        ]
        log_emission_matrix = [
            [math.log(p) for p in row] for row in self.emission_matrix
        ]

        # Initialize forward probability matrix
        alpha = [[-math.inf for _ in range(num_states)] for _ in range(num_obs)]

        # Initial probabilities
        for s in range(num_states):
            alpha[0][s] = log_initial_probs[s] + log_emission_matrix[s][observations[0]]

        # Compute probabilities for subsequent observations
        for t in range(1, num_obs):
            for j in range(num_states):
                alpha[t][j] = (
                    self.log_sum_exp(
                        [
                            alpha[t - 1][i] + log_transition_matrix[i][j]
                            for i in range(num_states)
                        ]
                    )
                    + log_emission_matrix[j][observations[t]]
                )

        # Total probability of the observation sequence
        log_prob = self.log_sum_exp(alpha[-1])

        return [[math.exp(a) for a in row] for row in alpha], math.exp(log_prob)

    def log_sum_exp(self, log_probs: List[float]) -> float:
        """Compute log(sum(exp(log_probs))) in a numerically stable way."""
        max_log_prob = max(log_probs)

        return max_log_prob + math.log(
            sum(math.exp(x - max_log_prob) for x in log_probs)
        )

    def backward(self, observations: List[int]) -> List[List[float]]:
        num_states = len(self.states)
        num_obs = len(observations)

        # Convert probabilities to log probabilities for numerical stability
        log_transition_matrix = [
            [math.log(p) for p in row] for row in self.transition_matrix
        ]

        log_emission_matrix = [
            [math.log(p) for p in row] for row in self.emission_matrix
        ]

        # Initialize backward matrix with log probabilities
        beta = [[-math.inf for _ in range(num_states)] for _ in range(num_obs)]

        for s in range(num_states):
            beta[num_obs - 1][s] = 0  # log(1)

        # Compute backward probabilities
        for t in range(num_obs - 2, -1, -1):
            for i in range(num_states):
                beta[t][i] = self.log_sum_exp(
                    [
                        beta[t + 1][j]
                        + log_transition_matrix[i][j]
                        + log_emission_matrix[j][observations[t + 1]]
                        for j in range(num_states)
                    ]
                )

        return [[math.exp(b) for b in row] for row in beta]

    def baum_welch_training(self, observations: List[int], n_iter=100):
        num_states = len(self.transition_matrix)
        num_obs = len(set(observations))

        for _ in range(n_iter):
            forward_probs, _ = self.forward(observations)
            backward_probs = self.backward(observations)

            # Initialize xi and gamma matrices
            xi = [
                [[0.0 for _ in range(len(observations) - 1)] for _ in range(num_states)]
                for _ in range(num_states)
            ]
            gamma = [[0.0 for _ in range(len(observations))] for _ in range(num_states)]

            for t in range(len(observations) - 1):
                denominator = sum(
                    forward_probs[t][i]
                    * self.transition_matrix[i][j]
                    * self.emission_matrix[j][observations[t + 1]]
                    * backward_probs[t + 1][j]
                    for i in range(num_states)
                    for j in range(num_states)
                )

                for i in range(num_states):
                    for j in range(num_states):
                        numerator = (
                            forward_probs[t][i]
                            * self.transition_matrix[i][j]
                            * self.emission_matrix[j][observations[t + 1]]
                            * backward_probs[t + 1][j]
                        )
                        xi[i][j][t] = numerator / denominator if denominator else 0
                        gamma[i][t] += xi[i][j][t]

            # Update transition matrix
            for i in range(num_states):
                for j in range(num_states):
                    numerator = sum(xi[i][j][t] for t in range(len(observations) - 1))
                    denominator = sum(gamma[i][t] for t in range(len(observations) - 1))
                    self.transition_matrix[i][j] = (
                        numerator / denominator if denominator else 0
                    )

            # Update emission matrix
            for j in range(num_states):
                for k in range(num_obs):
                    numerator = sum(
                        gamma[j][t]
                        for t in range(len(observations))
                        if observations[t] == k
                    )

                    denominator = sum(gamma[j][t] for t in range(len(observations)))
                    self.emission_matrix[j][k] = (
                        numerator / denominator if denominator else 0
                    )

    def gaussian_pdf(self, mean, std, x):
        """Calculate Gaussian probability density function without external libraries."""
        coefficient = 1 / (std * math.sqrt(2 * math.pi))
        exponent = -0.5 * ((x - mean) / std) ** 2
        return coefficient * math.exp(exponent)

    def continuous_emission_probability(self, mean, std, observation):
        """Compute emission probability for continuous distributions (Gaussian)."""
        return self.gaussian_pdf(mean, std, observation)

    def discrete_emission_probability(self, state, observation):
        """Compute emission probability for discrete distributions."""
        return self.emission_matrix[state][observation]


def test_hmm():
    # Instantiate the HMM with given parameters
    transition_matrix = [[0.5, 0.5], [0.4, 0.6]]  # Transition probabilities
    emission_matrix = [[0.2, 0.3, 0.3, 0.2], [0.3, 0.2, 0.2, 0.3]]  # Emission probabilities
    initial_state_probs = [0.5, 0.5]  # Initial state probabilities

    hmm = HMM(transition_matrix, emission_matrix, initial_state_probs)

    # Observations for testing
    observations_viterbi = [2, 2, 1, 0, 1, 3, 2, 0, 0]
    observations_forward_backward = [2, 2, 1, 0]

    # Test Viterbi Algorithm
    predicted_states = hmm.viterbi_algorithm(observations_viterbi)
    expected_viterbi_states = [0, 0, 0, 1, 1, 1, 1, 1, 1]
    if predicted_states == expected_viterbi_states:
        print("Viterbi Test passed!")
    else:
        print(
            "Viterbi Test failed. Predicted:",
            predicted_states,
            "Expected:",
            expected_viterbi_states,
        )

    forward_probs, forward_prob = hmm.forward(observations_forward_backward)
    expected_forward_prob = 0.0038432
    epsilon = 0.0001

    # Compare the calculated probability with the expected value
    if abs(forward_prob - expected_forward_prob) < epsilon:
        print("Forward algorithm test passed! Probability matches the expected value.")
    else:
        print(
            f"Forward algorithm test failed. Calculated probability: {forward_prob}, Expected probability: {expected_forward_prob}"
        )

    # Run the backward algorithm and compare with the expected result
    backward_probs_matrix = hmm.backward(observations_forward_backward)
    backward_prob = sum(
        hmm.initial_state_probs[i]
        * hmm.emission_matrix[i][observations_forward_backward[0]]
        * backward_probs_matrix[0][i]
        for i in range(len(hmm.states))
    )
    expected_backward_prob = 0.0038432
    epsilon = 0.0001

    # Compare the calculated probability with the expected value
    if abs(backward_prob - expected_backward_prob) < epsilon:
        print("Backward algorithm test passed! Probability matches the expected value.")
    else:
        print(
            f"Backward algorithm test failed. Calculated probability: {backward_prob}, Expected probability: {expected_backward_prob}"
        )


    # Instantiate the HMM with given parameters for Baum-Welch training
    transition_matrix = [[0.8, 0.1], [0.1, 0.8]]
    emission_matrix = [[0.7, 0.2, 0.1], [0.1, 0.2, 0.7]]
    initial_state_probs = [0.3, 0.7]

    hmm_bw = HMM(transition_matrix, emission_matrix, initial_state_probs)

    # Observations for testing Baum-Welch algorithm
    observations_bw = [
        1,
        2,
        2,
        1,
        2,
        1,
        2,
        1,
        1,
        2,
        0,
        2,
        2,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        2,
        0,
        1,
        0,
        0,
        0,
        1,
        2,
        2,
        1,
        2,
        1,
        1,
    ]

    # Expected parameters after Baum-Welch training
    expected_transition_matrix = [[0.8797, 0.1049], [0.0921, 0.8658]]
    expected_emission_matrix = [[0.6765, 0.2188, 0.1047], [0.0584, 0.4251, 0.5165]]


    # Function to compare matrices with a tolerance
    def matrices_are_close(matrix1, matrix2, tolerance=0.1):
        for row1, row2 in zip(matrix1, matrix2):
            for val1, val2 in zip(row1, row2):
                if abs(val1 - val2) > tolerance:
                    return False
        return True


    # Run Baum-Welch training and compare the results
    hmm_bw.baum_welch_training(observations_bw)

    # Compare the updated matrices with the expected values
    transition_matrix_test_passed = matrices_are_close(
        hmm_bw.transition_matrix, expected_transition_matrix
    )
    emission_matrix_test_passed = matrices_are_close(
        hmm_bw.emission_matrix, expected_emission_matrix
    )

    print(
        "Transition Matrix Test:", "Passed" if transition_matrix_test_passed else "Failed"
    )
    print("Updated Transition Matrix:", hmm_bw.transition_matrix)
    print("Emission Matrix Test:", "Passed" if emission_matrix_test_passed else "Failed")
    print("Updated Emission Matrix:", hmm_bw.emission_matrix)


    transition_matrix = [[0.5, 0.5], [0.4, 0.6]]  # Transition probabilities
    emission_matrix = [[0.2, 0.3, 0.3, 0.2], [0.3, 0.2, 0.2, 0.3]]  # Emission probabilities
    initial_state_probs = [0.5, 0.5]  # Initial state probabilities

    hmm = HMM(transition_matrix, emission_matrix, initial_state_probs)

    # Testing continuous emission probability
    mean_state_0, std_state_0 = 0.0, 1.0  # Example mean and standard deviation
    observation = 0.5  # Example observation
    continuous_prob = hmm.continuous_emission_probability(
        mean_state_0, std_state_0, observation
    )
    print("Continuous emission probability:", continuous_prob)

    # Testing discrete emission probability
    state, observation = 0, 2  # Example state and observation
    discrete_prob = hmm.discrete_emission_probability(state, observation)
    print("Discrete emission probability:", discrete_prob)


if __name__ == "__main__":
    test_hmm()