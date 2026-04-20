import json
import re
import os

class Lookup:
    @staticmethod
    def create_stack_json(data_file, json_file):
        """Parse the stack.dat file and create the 4D JSON."""
        with open(data_file, 'r') as f:
            lines = f.readlines()
        
        # Remove empty lines
        lines = [line for line in lines if line.strip()]
        
        # Find start of stacking energies
        start_idx = next(i for i, line in enumerate(lines) if 'STACKING ENERGIES' in line) + 1
        
        data_lines = lines[start_idx:]
        
        bases = ['A', 'C', 'G', 'U']
        stack_energy = {}
        
        for i in range(4):  # base_i
            block_start = i * 12
            data = data_lines[block_start + 8: block_start + 12]  # 4 data lines
            
            for row in range(4):  # base_jm1
                data_line = data[row]
                tokens = re.split(r'\s+', data_line.strip())
                if len(tokens) != 16:
                    continue  # skip if not 16 tokens
                for j in range(4):  # base_j
                    for col in range(4):  # base_ip1
                        val = tokens[j * 4 + col]
                        if val != '.':
                            base_i = bases[i]
                            base_j = bases[j]
                            base_ip1 = bases[col]
                            base_jm1 = bases[row]
                            stack_energy.setdefault(base_i, {}).setdefault(base_ip1, {}).setdefault(base_j, {}).setdefault(base_jm1, float(val))
        
        with open(json_file, 'w') as f:
            json.dump(stack_energy, f, indent=2)

    @staticmethod
    def create_loop_jsons(data_file, hairpin_file, bulge_file, internal_file):
        """Parse the loop.dat file and create the 1D JSONs for hairpin, bulge, internal."""
        with open(data_file, 'r') as f:
            lines = f.readlines()
        
        hairpin = {}
        bulge = {}
        internal = {}
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('SIZE') or line.startswith('---') or line.startswith('DEST'):
                continue
            tokens = re.split(r'\s+', line)
            if len(tokens) == 4:
                size, int_energy, bulge_energy, hairpin_energy = tokens
                size = str(int(size))
                if int_energy != '.':
                    internal[size] = float(int_energy)
                if bulge_energy != '.':
                    bulge[size] = float(bulge_energy)
                if hairpin_energy != '.':
                    hairpin[size] = float(hairpin_energy)
        
        with open(hairpin_file, 'w') as f:
            json.dump(hairpin, f, indent=2)
        with open(bulge_file, 'w') as f:
            json.dump(bulge, f, indent=2)
        with open(internal_file, 'w') as f:
            json.dump(internal, f, indent=2)

    def __init__(self, json_path=None):
        # Default stack JSON path inside the algos package
        if json_path is None:
            json_path = os.path.join(os.path.dirname(__file__), 'stack.json')

        data_file = os.path.join(os.path.dirname(__file__), '..', 'datasets', 'stack.dat')
        if not os.path.exists(json_path):
            self.create_stack_json(data_file, json_path)

        # Load stack energy
        with open(json_path, 'r') as f:
            self.stack_energy = json.load(f)

        # Create loop JSONs if not exist
        loop_data_file = os.path.join(os.path.dirname(__file__), '..', 'datasets', 'loop.dat')
        dir_path = os.path.dirname(json_path)
        hairpin_file = os.path.join(dir_path, 'hairpin.json')
        bulge_file = os.path.join(dir_path, 'bulge.json')
        internal_file = os.path.join(dir_path, 'internal.json')

        if not os.path.exists(hairpin_file) or not os.path.exists(bulge_file) or not os.path.exists(internal_file):
            self.create_loop_jsons(loop_data_file, hairpin_file, bulge_file, internal_file)

        # Load loop energies
        with open(hairpin_file, 'r') as f:
            self.hairpin_energy = json.load(f)
        with open(bulge_file, 'r') as f:
            self.bulge_energy = json.load(f)
        with open(internal_file, 'r') as f:
            self.internal_energy = json.load(f)

    def stack(self, base_i, base_ip1, base_j, base_jm1):
        """Get the stacking energy for the given bases.
        
        Args:
            base_i: base at position i
            base_ip1: base at position i+1
            base_j: base at position j
            base_jm1: base at position j-1
            
        Returns:
            float: stacking free energy, 0.0 if not found
        """
        return self.stack_energy.get(base_i, {}).get(base_ip1, {}).get(base_j, {}).get(base_jm1, 0.0)

    def hairpin(self, size):
        """Get the hairpin loop energy for the given size."""
        return self.hairpin_energy.get(str(size), 0.0)

    def bulge(self, size):
        """Get the bulge loop energy for the given size."""
        return self.bulge_energy.get(str(size), 0.0)

    def internal(self, size):
        """Get the internal loop energy for the given size."""
        return self.internal_energy.get(str(size), 0.0)

    def multiloop(self, p, u):
        """Get the multiloop energy for p inner pairs and u unpaired bases.

        Energy = a + b*p + c*u where a=2, b=42, c=0.
        """
        a = 2
        b = 42
        c = 0
        return a + b * p + c * u

if __name__ == '__main__':
    lookup_instance = Lookup()
    print('stack(A,U,U,G)=', lookup_instance.stack('A', 'U', 'U', 'G'))
    print('hairpin(3)=', lookup_instance.hairpin(3))
    print('bulge(3)=', lookup_instance.bulge(3))
    print('internal(4)=', lookup_instance.internal(4))
    print('multiloop(4)=', lookup_instance.multiloop(4,2))
