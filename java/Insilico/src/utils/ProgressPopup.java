package utils;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import utils.ProgressBar;

public class ProgressPopup {
	JFrame frame;
	ProgressBar pb;
	public ProgressPopup(boolean isUnix, int N){
		if(!isUnix){
			    this.frame = new JFrame("ProgressBar");
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			    //Create and set up the content pane.
			    this.pb = new ProgressBar(N);
			    JComponent newContentPane = pb;
			    newContentPane.setOpaque(true); //content panes must be opaque
			    frame.setContentPane(newContentPane);
			    //Display the window.
			    frame.pack();
			    frame.setVisible(true);
		}
	}
	
	public void update(int k){
		this.pb.change(k);
	}
	
	public void message(int k, String message){
		this.pb.print(k, message);
	}
	
	public void close(){
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                frame.setVisible(false);
                frame.dispose();
            }
        });
	}
}
