package GUI;
import static GUI.Utils.*;

import java.awt.BorderLayout;
import java.awt.Dimension;

import javax.swing.JPanel;
import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;

public final class TreePanel extends Panel {
	
	private JTree tree;
	
	TreePanel(){
		super("Arborescence des fichiers", yellow);
		setBackground(darkGrey);
		initTree();
	}

	protected void initTree() {
		DefaultMutableTreeNode root=new DefaultMutableTreeNode("Genbank");  
	    DefaultMutableTreeNode archaea=new DefaultMutableTreeNode("Archaea"); 
	    DefaultMutableTreeNode bacteria=new DefaultMutableTreeNode("Bacteria");  
	    DefaultMutableTreeNode group=new DefaultMutableTreeNode("Group");  
	    DefaultMutableTreeNode group2=new DefaultMutableTreeNode("Group");
	    DefaultMutableTreeNode subgroup=new DefaultMutableTreeNode("SubGroup"); 
	    DefaultMutableTreeNode organism=new DefaultMutableTreeNode("Organism"); 
	    root.add(archaea);  
	    root.add(bacteria);  
	    archaea.add(group); 
	    bacteria.add(group2);
	    group.add(subgroup);
	    subgroup.add(organism);
	    DefaultMutableTreeNode red=new DefaultMutableTreeNode("Organism 1");  
	    DefaultMutableTreeNode blue=new DefaultMutableTreeNode("Organism 2");  
	    organism.add(red); organism.add(blue);  
	    tree = new JTree(root);  
        tree.setShowsRootHandles(true);
	    add(tree, BorderLayout.CENTER);
	    tree.setPreferredSize(new Dimension(500, 600));
	    tree.setBackground(darkGrey);
	    
	    DefaultTreeCellRenderer renderer = (DefaultTreeCellRenderer) tree.getCellRenderer();
        renderer.setTextSelectionColor(white);
        renderer.setBackgroundSelectionColor(grey);
        renderer.setBackgroundNonSelectionColor(darkGrey);
        renderer.setTextNonSelectionColor(white);
        renderer.setBorderSelectionColor(black);
        renderer.setLeafIcon(null);
	}
	
	protected void createComponents() {
		// TODO Auto-generated method stub
		
	}

	protected void addComponents() {
		// TODO Auto-generated method stub
		
	}

	protected void styleComponents() {
		
	}
}
