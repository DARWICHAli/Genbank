package GUI;

import static GUI.Utils.*;

import java.awt.*;

import javax.swing.*;
import javax.swing.border.Border;

public abstract class Panel extends JPanel {
	
	Panel(String title, Color color) {
		final TitlePanel titlePanel = new TitlePanel(title, 22, new Color(0, 120, 130), white);
		//titlePanel.setBorder(BorderFactory.createLineBorder(yellow));
		this.setLayout(new BorderLayout());
        this.add(titlePanel, BorderLayout.NORTH);
        titlePanel.setPreferredSize(new Dimension(400, 50));
        createComponents();
        addComponents();
        styleComponents();
	}
	
	
	Panel(String title) {
		final TitlePanel titlePanel = new TitlePanel(title, 22, new Color(50, 60, 70), white);
		this.setLayout(new BorderLayout());
        this.add(titlePanel, BorderLayout.NORTH);
        this.setBackground(lightGrey);
        createComponents();
        addComponents();
        styleComponents();
	}
	
	Panel() {
		this.setLayout(new BorderLayout());
        this.setBackground(lightGrey);
        createComponents();
        addComponents();
        styleComponents();
	}
	
	protected abstract void createComponents();
	protected abstract void addComponents();
	protected abstract void styleComponents();
}
