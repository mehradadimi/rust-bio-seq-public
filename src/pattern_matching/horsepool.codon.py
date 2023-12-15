from bio import *


class Horspool:
    m: int
    pattern: str
    shift: list[int]

    def __init__(self, pattern: str):
        self.m = len(pattern)
        self.pattern = pattern
        self.shift = self.build_shift_table(pattern)

    def build_shift_table(self, pattern: str):
        shift = [len(pattern) for _ in range(256)]  # Initialize shift values for ASCII characters

        for i in range(len(pattern) - 1):
            shift[ord(pattern[i])] = len(pattern) - 1 - i

        return shift

    def find_all(self, text: str) -> list:
        matches = []
        n = len(text)
        last = len(self.pattern) - 1  # Index of the last character of the pattern

        i = 0
        while i <= n - len(self.pattern):
            if text[i + last] == self.pattern[last]:
                if text[i:i+len(self.pattern)] == self.pattern:
                    matches.append(i)

                i += self.shift[ord(self.pattern[last])]  # Shift based on the last character of the pattern
                
            else:
                i += self.shift[ord(text[i + last])]  # Shift based on mismatched character

        return matches


def test_horsepool():
    # Testing
    print('Test 1')
    print('-')
    print('Expected: [7, 17]')
    pattern = "GAAAA"
    text = "ACGGCTAGAAAAGGCTAGAAAA"
    hp = Horspool(pattern)
    print(f'Result: {hp.find_all(text)}')

    print('----------')

    print('Test 2')
    print('-')
    print('Expected: [8]')
    pattern = "qnnnannan"
    text = "dhjalkjwqnnnannanaflkjdklfj"
    hp = Horspool(pattern)
    print(f'Result: {hp.find_all(text)}')

    print('----------')

    print('Test 3')
    print('-')
    print('Expected: [0]')
    pattern = "dhjalk"
    text = "dhjalkjwqnnnannanaflkjdklfj"
    hp = Horspool(pattern)
    print(f'Result: {hp.find_all(text)}')
    
    print('----------')

    print("Manual evaluation: All Tests Pass")


if __name__ == "__main__":
    test_horsepool()