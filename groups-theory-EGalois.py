import itertools

class Groupe:
    def __init__(self, elements, operation):
        self.elements = elements
        self.operation = operation
        self.identity = self.find_identity()
        self.inverses = self.find_inverses()

    def find_identity(self):
        for e in self.elements:
            if all(self.operation(e, a) == a and self.operation(a, e) == a for a in self.elements):
                return e
        return None

    def find_inverses(self):
        inverses = {}
        for a in self.elements:
            for b in self.elements:
                if self.operation(a, b) == self.identity and self.operation(b, a) == self.identity:
                    inverses[a] = b
                    break
        return inverses

    def is_group(self):
        # Check closure
        for a in self.elements:
            for b in self.elements:
                if self.operation(a, b) not in self.elements:
                    return False
        
        # Check associativity
        for a in self.elements:
            for b in self.elements:
                for c in self.elements:
                    if self.operation(self.operation(a, b), c) != self.operation(a, self.operation(b, c)):
                        return False
        
        # Check identity element
        if self.identity is None:
            return False
        
        # Check inverses
        if len(self.inverses) != len(self.elements):
            return False
        
        return True

def permutation_multiplication(p1, p2):
    return tuple(p1[i-1] for i in p2)

# Exemple de groupe de permutations
elements = list(itertools.permutations([1, 2, 3]))
operation = permutation_multiplication

groupe_permutations = Groupe(elements, operation)

print("Est-ce un groupe ? :", groupe_permutations.is_group())
print("Élément neutre :", groupe_permutations.identity)
print("Inverses :", groupe_permutations.inverses)

# Illustration graphique des permutations
import matplotlib.pyplot as plt
import networkx as nx

def draw_permutation_graph(elements, operation):
    G = nx.DiGraph()
    for e in elements:
        G.add_node(e)
    for e1 in elements:
        for e2 in elements:
            G.add_edge(e1, e2, label=str(operation(e1, e2)))
    
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_size=2000, node_color='lightblue', font_size=10)
    edge_labels = nx.get_edge_attributes(G, 'label')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red')
    plt.title("Groupe de Permutations")
    plt.show()

draw_permutation_graph(elements, operation)