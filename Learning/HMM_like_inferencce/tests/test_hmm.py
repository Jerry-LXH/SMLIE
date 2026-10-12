import unittest
import numpy as np
from src.hmm import forward_backward, brute_force, poisson_log_emissions
class TestHMM(unittest.TestCase):
    def test_exact_and_marginals(self):
        pi=np.array([.7,.3]); A=np.array([[.85,.15],[.2,.8]])
        b=poisson_log_emissions([0,3,1,5],[.5,4])
        ll,g,x,f=forward_backward(b,pi,A); ref,gr,*_=brute_force(b,pi,A)
        self.assertAlmostEqual(ll,ref,places=12)
        np.testing.assert_allclose(g,gr,atol=1e-12)
        np.testing.assert_allclose(x.sum(axis=2),g[:-1],atol=1e-12)
        np.testing.assert_allclose(x.sum(axis=1),g[1:],atol=1e-12)
        np.testing.assert_allclose(f.sum(axis=1),1)
    def test_long_sequence(self):
        ll,g,*_=forward_backward(poisson_log_emissions([1]*2000,[.5,4]),np.array([1.,0.]),np.array([[.9,.1],[0.,1.]]))
        self.assertTrue(np.isfinite(ll)); np.testing.assert_allclose(g.sum(axis=1),1,atol=1e-9)
    def test_one_frame(self):
        _,g,x,_=forward_backward(poisson_log_emissions([2],[1.,3.]),np.array([.5,.5]),np.eye(2))
        self.assertEqual(x.shape,(0,2,2)); self.assertAlmostEqual(g.sum(),1)
if __name__=='__main__': unittest.main()
