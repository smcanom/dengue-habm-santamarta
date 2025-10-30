#!/usr/bin/env python3
"""
test_params.py

Quick verification that parameter parsing and perturbation logic works correctly.
"""

import re
from pathlib import Path


def rewrite_props(props_path: Path, updates: dict) -> None:
    """Rewrite a properties file with updated key=value pairs."""
    text = props_path.read_text(encoding='utf-8', errors='ignore')
    for key, val in updates.items():
        pattern = re.compile(rf"^\s*{re.escape(key)}\s*=.*$", flags=re.MULTILINE)
        if pattern.search(text):
            text = pattern.sub(f"{key}={val}", text)
        else:
            if not text.endswith("\n"):
                text += "\n"
            text += f"{key}={val}\n"
    props_path.write_text(text, encoding='utf-8')


def read_param(props_path: Path, key: str) -> str:
    """Read a parameter value from props file."""
    text = props_path.read_text(encoding='utf-8', errors='ignore')
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s*(.*)$", flags=re.MULTILINE)
    match = pattern.search(text)
    if match:
        return match.group(1).strip()
    return ""


def test_parameter_perturbation():
    """Test that parameter perturbation works correctly."""
    
    print("=" * 60)
    print("Parameter Perturbation Test")
    print("=" * 60)
    
    # Create a test props file
    test_props = Path("test_model.props")
    
    # Write baseline
    test_props.write_text("""
# Test properties
sigma_M=0.3
sigma_H=3.0
z=0.3
r=0.6
C=30.0
beta_mh=0.10
beta_hm=0.10
config_id=baseline
perturb_param=baseline
perturb_delta=0.0
replicate_id=0
base_seed=12345
""", encoding='utf-8')
    
    print("\n1. Baseline Configuration:")
    print(f"   sigma_M = {read_param(test_props, 'sigma_M')}")
    print(f"   z = {read_param(test_props, 'z')}")
    print(f"   config_id = {read_param(test_props, 'config_id')}")
    
    # Test perturbation
    print("\n2. Applying +10% perturbation to z:")
    
    # Read original
    z_orig = float(read_param(test_props, 'z'))
    
    # Apply perturbation (simulating what Python driver does)
    z_pert = z_orig * 1.10
    
    updates = {
        'z': str(z_pert),
        'config_id': 'z_p10',
        'perturb_param': 'z',
        'perturb_delta': '0.10',
        'replicate_id': '1'
    }
    
    rewrite_props(test_props, updates)
    
    print(f"   z_original = {z_orig}")
    print(f"   z_perturbed = {read_param(test_props, 'z')}")
    print(f"   Expected = {z_pert}")
    print(f"   config_id = {read_param(test_props, 'config_id')}")
    print(f"   perturb_param = {read_param(test_props, 'perturb_param')}")
    print(f"   replicate_id = {read_param(test_props, 'replicate_id')}")
    
    # Verify
    assert abs(float(read_param(test_props, 'z')) - z_pert) < 1e-9, "Perturbation failed!"
    assert read_param(test_props, 'config_id') == 'z_p10', "Config ID not updated!"
    
    print("\n3. Testing CRN (Common Random Numbers):")
    base_seed = int(read_param(test_props, 'base_seed'))
    replicate_id = int(read_param(test_props, 'replicate_id'))
    run_seed = base_seed + replicate_id
    
    print(f"   base_seed = {base_seed}")
    print(f"   replicate_id = {replicate_id}")
    print(f"   run_seed = base_seed + replicate_id = {run_seed}")
    print(f"   (For CRN: same replicate_id across all configs → same random stream)")
    
    # Clean up
    test_props.unlink()
    
    print("\n" + "=" * 60)
    print("✓ All parameter tests passed!")
    print("=" * 60)
    print("\nNext Steps:")
    print("  1. Compile the model: make clean && make")
    print("  2. Run sensitivity analysis: python run_sensitivity.py --mode run --replicates 3")
    print("  3. Aggregate results: python run_sensitivity.py --mode aggregate")
    print("=" * 60)


if __name__ == '__main__':
    test_parameter_perturbation()

