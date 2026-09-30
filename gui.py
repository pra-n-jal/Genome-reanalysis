import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import subprocess
import os

class PipelineApp:
    def __init__(self, root):
        self.root = root
        # TODO: maybe add an icon later?
        self.root.title("My Pipeline Runner v1.2 (DO NOT CLOSE DURING RUN)")
        self.root.geometry("750x650")
        
        # --- UI LAYOUT STUFF ---
        
        # 1. Input File Frame
        step1_frame = tk.LabelFrame(root, text="Step 1: Pick the data file", padx=10, pady=10)
        step1_frame.pack(fill="x", padx=10, pady=5)
        
        self.file_path_var = tk.StringVar()
        # width 80 is big enough for most windows paths
        tk.Entry(step1_frame, textvariable=self.file_path_var, width=80).pack(side="left", padx=5)
        tk.Button(step1_frame, text="Browse...", command=self.browse_file).pack(side="left")
        
        # 2. Script Selection Frame
        step2_frame = tk.LabelFrame(root, text="Step 2: Choose what to do", padx=10, pady=10)
        step2_frame.pack(fill="x", padx=10, pady=5)
        
        # my scripts (format: nice name, actual file, command line flag for input)
        self.my_scripts = [
            ("Signature: Build C4A", "src/signature/build_c_4_a_complement_signature.py", "--input"),
            ("Signature: Convert to Entrez", "src/signature/convert_signature_to_entr.py", "--input"),
            ("LINCS: Extract Metadata", "src/lincs/extract_lincs_metadata.py", "--gctx"),
            ("LINCS: Compute Connectivity", "src/lincs/compute_connectivity_lincs.py", "--sig"),
            ("Chemprop: Filter and Cluster", "src/chemprop/filter_and_cluster.py", "--input")
        ]
        
        self.script_var = tk.StringVar(value=self.my_scripts[0][0])
        self.dropdown = ttk.Combobox(
            step2_frame, 
            textvariable=self.script_var, 
            values=[s[0] for s in self.my_scripts], 
            width=60, 
            state="readonly"
        )
        self.dropdown.pack(side="left", padx=5)
        
        # 3. Extra Arguments Frame (mostly for --out stuff)
        step3_frame = tk.LabelFrame(root, text="Step 3: Extra Arguments (like --out something)", padx=10, pady=10)
        step3_frame.pack(fill="x", padx=10, pady=5)
        
        self.args_var = tk.StringVar()
        tk.Entry(step3_frame, textvariable=self.args_var, width=80).pack(side="left", padx=5)
        
        # --- Action Buttons ---
        btn_frame = tk.Frame(root)
        btn_frame.pack(fill="x", padx=10, pady=10)
        
        # The main run button
        tk.Button(
            btn_frame, text="RUN PIPELINE SCRIPT", command=self.run_script, 
            bg="#d0d0d0", font=("Arial", 12, "bold"), pady=10, width=40
        ).pack(side="left", padx=5)
        
        # quick test button so i don't have to type it out every time
        tk.Button(
            btn_frame, text="LOAD TEST EXAMPLE", command=self.load_test_example, 
            bg="#a0c0ff", font=("Arial", 10, "bold"), pady=10
        ).pack(side="right", padx=5)
        
        # 4. Console Output so we can see print() statements
        console_frame = tk.LabelFrame(root, text="Terminal Output", padx=10, pady=10)
        console_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        # gotta add a scrollbar or we can't see the top errors
        scrollbar = tk.Scrollbar(console_frame)
        scrollbar.pack(side="right", fill="y")
        
        # green on black looks cooler
        self.console = tk.Text(console_frame, bg="black", fg="#00ff00", font=("Consolas", 10), yscrollcommand=scrollbar.set)
        self.console.pack(fill="both", expand=True)
        scrollbar.config(command=self.console.yview)
        
        print("GUI started up okay.")
        
    def browse_file(self):
        # find the root dir of this script so we default to data/
        my_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(my_dir, "data")
        
        # fallback just in case data/ got deleted
        if not os.path.exists(data_dir):
            data_dir = my_dir
            
        filename = filedialog.askopenfilename(initialdir=data_dir)
        if filename:
            self.file_path_var.set(filename)
            
    def run_script(self):
        selected_name = self.script_var.get()
        infile = self.file_path_var.get()
        
        target_script = ""
        target_flag = ""
        
        # loop to find the script path based on the nice name
        for nice_name, actual_path, flag in self.my_scripts:
            if nice_name == selected_name:
                target_script = actual_path
                target_flag = flag
                break
                
        if infile == "":
            messagebox.showwarning("Hold on!", "You forgot to select an input file first!")
            return
            
        # hack to make sure it runs even if we launched from another folder
        my_dir = os.path.dirname(os.path.abspath(__file__))
        absolute_script = os.path.join(my_dir, target_script)

        # Build up the python command
        cmd = ["python", absolute_script, target_flag, infile]
        
        # add whatever the user typed in the extra args box
        extras = self.args_var.get().strip()
        if extras != "":
            # split by space to pass them properly to subprocess
            cmd.extend(extras.split())
            
        # dump the command to the console so we know what's running
        self.console.insert(tk.END, "\n========================================\n")
        self.console.insert(tk.END, f"> RUNNING: {' '.join(cmd)}\n")
        self.console.see(tk.END)
        self.root.update() # force tkinter to draw the text before freezing
        
        # try running it
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.stdout:
                self.console.insert(tk.END, res.stdout)
            if res.stderr:
                self.console.insert(tk.END, "\n[LOGS/ERRORS]\n" + res.stderr)
        except Exception as e:
            self.console.insert(tk.END, f"\nWhoops, it crashed: {str(e)}\n")
            
        self.console.insert(tk.END, "========================================\n")
        self.console.see(tk.END)

    def load_test_example(self):
        # fill in the blanks with a good file to make testing faster
        my_dir = os.path.dirname(os.path.abspath(__file__))
        test_f = os.path.join(my_dir, "data", "intermediate", "C4A_complement_signature.tsv")
        
        self.file_path_var.set(test_f)
        
        if os.path.exists(test_f):
            self.script_var.set("Signature: Convert to Entrez")
            # put the output right into outputs/ so it doesn't clutter
            out_path = os.path.join(my_dir, 'data', 'outputs', 'test_entrez_output.tsv')
            self.args_var.set(f"--out {out_path}")
            
            self.console.insert(tk.END, "\n[DEBUG] Loaded the test setup! Hit the RUN button to test it out.\n")
        else:
            self.console.insert(tk.END, f"\n[ERROR] Wait, where is {test_f}??\n")
        
        self.console.see(tk.END)

if __name__ == "__main__":
    main_window = tk.Tk()
    app = PipelineApp(main_window)
    main_window.mainloop()
