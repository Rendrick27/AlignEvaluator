import os
import toytree
import toyplot
import toyplot.svg


def TreeGenerator(infile, outfile):
    """
    Generate a phylogenetic tree and save it as an SVG file.

    Args:
        infile (str): Path to the input file containing the tree.
        outfile (str): Path to save the output SVG file.
    """
    # Check if input file exists
    if not os.path.isfile(infile):
        raise FileNotFoundError(f"Input file '{infile}' not found.")

    # Read the Newick tree string from the input file
    with open(infile, 'r') as f:
        newick = f.read()

    # Create a Toytree object from the Newick string
    tree = toytree.tree(newick)

    # Remove nodes with only one child
    tree = remove_single_child_nodes(tree)

    # Get the number of taxa in the tree
    num_taxa = len(tree.get_tip_labels())

    # Calculate width and height based on the number of taxa
    width = max(2000, num_taxa * 50)  # Adjust the minimum width as needed
    height = max(1000, num_taxa * 30)  # Adjust the minimum height as needed

    # Get the support values for internal nodes
    support_values = tree.get_node_values('support')

    # Round and format the support values
    labels = [f"{round(float(s), 2)}" if s else "" for s in support_values]

    # Create the phylogenetic tree plot with support labels
    canvas, axes, mark = tree.draw(node_labels=labels, node_sizes=30,
                                   width=width, height=height)

    # Save the phylogenetic tree plot as an SVG file
    toyplot.svg.render(canvas, outfile)


def remove_single_child_nodes(tree):
    """
    Remove nodes with only one child from the tree.
    """
    while True:
        singles = [node.idx for node in tree.idx_dict.values()
                   if len(node.children) == 1]
        if not singles:
            break
        tree = tree.drop(treenodes=singles, axis=0)
    return tree
