import unittest

class TestSomething(unittest.TestCase):

	def test_some_method(self):
		self.assertEqual(2 + 2, 3)

if __name__ == "__main__":
    unittest.main()
