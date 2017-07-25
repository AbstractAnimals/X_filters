#!/usr/bin/python

import X_filtering as Xf
import os

""""
Script to test the basic functionality of functions in the X_filtering.py script.
"""
scriptdir = os.path.dirname(__file__)
print('Script dir %s: ' % scriptdir)

def test_findgenders():
    """
    Test the find genders function
    """
    x = ['AMf1', 'AMm1', 'AMf2']
    offset = 2
    males, females = Xf.find_genders(x, offset)
    assert(males == [3])
    assert(females == [2, 4])

    r_males, r_females = Xf.find_genders(x, offset, reverse=True)
    assert(r_females == [3])
    assert(r_males == [2, 4])

def test_df_totals():
    totals = Xf.total_read_dp_per_individual(scriptdir + '/test.vcf', 9, 20)
    assert(totals[9] == 233)
    assert(totals[10] == 450)
    assert(totals[11] == 174)

if __name__ == '__main__':
    test_df_totals()

