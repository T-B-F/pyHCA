from __future__ import division, absolute_import, print_function

import os, sys, shutil
import numpy as np
import unittest

from pyHCA import HCA


class TestSegHCA(unittest.TestCase):
    
    def setup(self):
        self.check_sequence = "NGRHTGFGRTCCDKGADHLKGEGHCCITLAKRGYFPCEPWCTLLFALNMFNMQNMMRQQFSDDHNNMGRLCQQTTHRFPFNSDNKEEYIWLYKVQRLGAW"
        self.check_domains = {0: [26, 70, 0.0021107348166291562, -0.11363636363636363],
                              1: [87, 100, 2.5160743828633869e-05, 1.6153846153846154]}
        self.check_clusters = {0: [ 6,  7,  np.array([1])],
                               1: [18, 19,  np.array([1])],
                               2: [26, 29,  np.array([1, 0, 1])],
                               3: [33, 35,  np.array([1, 1])],
                               4: [39, 60,  np.array([1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1])],
                               5: [66, 70,  np.array([1, 0, 0, 1])],
                               6: [77, 78,  np.array([1])],
                               7: [79, 80,  np.array([1])],
                               8: [87, 100, np.array([1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1])],
                               }
    def assert_domain(self, i, dom):
        self.assertEqual(self.check_domains[i][0], dom.start)
        self.assertEqual(self.check_domains[i][1], dom.stop)
        assert(abs(self.check_domains[i][2] - dom.pvalue) < 1e-6)
        assert(abs(self.check_domains[i][3] - dom.score) < 1e-6)
        
    def test_domains(self):
        self.setup()
        hca = HCA(seq=self.check_sequence)
        hca.segments()
        domains = hca.get_domains()
        self.assertEqual(len(domains), len(self.check_domains))
        for i, dom in enumerate(domains):
            self.assert_domain(i, dom)
    
    def assert_cluster(self, i, cluster):
        self.assertEqual(self.check_clusters[i][0], cluster.start)
        self.assertEqual(self.check_clusters[i][1], cluster.stop)
        assert(np.all(self.check_clusters[i][2] == cluster.hydro_cluster))
        
    def test_clusters(self):
        self.setup()
        hca = HCA(seq=self.check_sequence)
        hca.segments()
        clusters = hca.get_clusters()
        self.assertEqual(len(clusters), len(self.check_clusters))
        for i, cluster in enumerate(clusters):
            self.assert_cluster(i, cluster)
    
    
if __name__ == "__main__":
    unittest.main()