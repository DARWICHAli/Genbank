package GUI;
import static GUI.Utils.*;

import java.awt.BorderLayout;

import javax.swing.*;

public final class InfoPanel extends Panel {
	
	private JTextPane textPane; 
	
	InfoPanel(){
		super("Information");
		
		createComponents();
		addComponents();
	}

	@Override
	protected void createComponents() {
		 textPane = new JTextPane();
		 textPane.setSize(WIDTH, HEIGHT*8/10);
		 textPane.setText("info section here, Maybe add some image slides for entertainment:-)");
		 textPane.setBackground(grey);
		 textPane.setForeground(white);
	}

	@Override
	protected void addComponents() {
		add(textPane, BorderLayout.CENTER);
	}

	@Override
	protected void styleComponents() {
		
	}
}
