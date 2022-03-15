import unittest
import numpy as np
import math

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

class MyTestCase(unittest.TestCase):
    def test_something(self):
       time = [1,2,3]
       a = find_nearest(time,2)
       print(a)
       time[a-1]

if __name__ == '__main__':
    unittest.main()
