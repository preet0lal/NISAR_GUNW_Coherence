"""Unit tests for NISAR GUNW native swath and exception masks.

Run from repo root:
    python -m unittest discover -s tests -v
"""
import unittest
import numpy as np
from nisar_gunw_coherence_ease2 import (
    _gunw_mask_components, _gunw_mask_is_valid, _inward_buffer,
)


class MaskDecoding(unittest.TestCase):
    def test_legacy_wrs_preserves_water_and_low_coherence_support(self):
        m = np.array([[11, 111, 10, 1, 255]], dtype=np.uint8)
        self.assertEqual(_gunw_mask_is_valid(m).tolist(),
                         [[True, True, False, False, False]])

    def test_packed_exception_bytes(self):
        m = np.array([[11, (1 << 8) | 11, (2 << 16) | 11,
                       (3 << 16) | (1 << 8) | 11, 0, 255]],
                     dtype=np.uint32)
        swath, exceptions = _gunw_mask_components(m)
        self.assertEqual(swath.tolist(),
                         [[True, True, True, True, False, False]])
        self.assertEqual(exceptions.tolist(),
                         [[True, False, False, False, True, True]])
        self.assertEqual(_gunw_mask_is_valid(m).tolist(),
                         [[True, False, False, False, False, False]])
        self.assertEqual(_gunw_mask_is_valid(m, use_exception_mask=False).tolist(),
                         [[True, True, True, True, False, False]])

    def test_negative_fill_sign_preserving(self):
        m = np.array([-1, 11], dtype=np.int32)
        self.assertEqual(_gunw_mask_is_valid(m).tolist(), [False, True])

    def test_edge_buffer(self):
        m = np.ones((7,7), dtype=bool)
        m[3,3] = False
        out = _inward_buffer(m, radius=2)
        self.assertFalse(out[3,3])
        self.assertFalse(out[3,5])
        self.assertTrue(out[0,0] == False)  # image outer boundary is eroded


if __name__ == '__main__':
    unittest.main()
