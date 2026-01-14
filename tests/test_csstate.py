from sandlercorrespondingstates.csstate import CSState
import pytest
from unittest import TestCase

class TestCSState(TestCase):

    def test_properties_full_supercrit(self):
        state = CSState(
            T=700.0,  # K
            P=5,       # MPa
            Tc=600.0,           # K
            Pc=4.0,             # MPa
            round=3
        )
        assert abs(state.Tr - (700.0 / 600.0)) < 1e-6
        assert abs(state.Pr - (5e6 / (4.0e6))) < 1e-6
            
        assert state.Z is not None
        self.assertEqual(state.Z.shape, (1,))
        assert state.h_dep is not None
        self.assertEqual(state.h_dep.shape, (1,))
        assert state.s_dep is not None
        self.assertEqual(state.s_dep.shape, (1,))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            CSState(
                Tc=600.0,
                Pc=4.0
            )
        
        with self.assertRaises(ValueError):
            CSState(
                T=700.0,
                P=5.0,
                x=0.5,
                Tc=600.0,
                Pc=4.0
            )
        
        with self.assertRaises(ValueError):
            CSState(
                T=700.0,
                x=0.5,
                Pc=4.0,
                Tc=600.0
            )
        
        with self.assertRaises(ValueError):
            CSState(
                P=5.0,
                x=0.5,
                Tc=600.0,
                Pc=4.0
            )

    def test_properties_subcrit_T_supercrit_P(self):
        state = CSState(
            T=500.0,  # K
            P=4.5,       # MPa
            Tc=600.0,           # K
            Pc=4.0,             # MPa
            round=3
        )
        assert abs(state.Tr - (500.0 / 600.0)) < 1e-6
        assert abs(state.Pr - (4.5e6 / (4.0e6))) < 1e-6
            
        assert state.Z is not None
        self.assertEqual(state.Z.shape, (1,))
        assert state.h_dep is not None
        self.assertEqual(state.h_dep.shape, (1,))
        assert state.s_dep is not None
        self.assertEqual(state.s_dep.shape, (1,))

    def test_properties_subcrit_P_supercrit_T(self):
        state = CSState(
            T=700.0,  # K
            P=3.5,       # MPa
            Tc=600.0,           # K
            Pc=4.0,             # MPa
            round=3
        )
        assert abs(state.Tr - (700.0 / 600.0)) < 1e-6
        assert abs(state.Pr - (3.5e6 / (4.0e6))) < 1e-6
            
        assert state.Z is not None
        self.assertEqual(state.Z.shape, (1,))
        assert state.h_dep is not None
        self.assertEqual(state.h_dep.shape, (1,))
        assert state.s_dep is not None
        self.assertEqual(state.s_dep.shape, (1,))

    def test_properties_two_phase(self):
        state = CSState(
            T=520.0,  # K
            x = 1.0,       # vapor phase
            Tc=600.0,           # K
            Pc=4.0,             # MPa
            round=3
        )
        assert abs(state.Tr - (520.0 / 600.0)) < 1e-6
            
        assert state.Z is not None
        self.assertEqual(state.Z.shape, (2,))
        assert state.h_dep is not None
        self.assertEqual(state.h_dep.shape, (2,))
        assert state.s_dep is not None
        self.assertEqual(state.s_dep.shape, (2,))
    
    def test_husv(self):
        state = CSState(
            T=400.0,  # K
            P=5,       # MPa
            round=3
        ).set_compound('methane')
        h, u, s, v = state.h, state.u, state.s, state.v
        self.assertEqual(h.shape, (1,))
        self.assertEqual(u.shape, (1,))
        self.assertEqual(s.shape, (1,))
        self.assertEqual(v.shape, (1,))