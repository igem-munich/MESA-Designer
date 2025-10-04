# MESA Design Tool

### 

### **1\. User Interface (UI)**

The user interface would be the core of the tool, guiding the user through a step-by-step process of selecting each MESA component. It would be highly visual and interactive.

**a. Step 1: Select the Target Ligand** The user would first choose the extracellular molecule they want to detect. The interface could provide a searchable list of common biomarkers, cytokines, or other molecules. This selection is crucial as it determines which ligand-binding domains are relevant.

**b. Step 2: Build the Receptor Chains** The user would then build the two MESA receptor chains, often referred to as the "Protease Chain" (PC) and the "Target Chain" (TC). For each chain, the user would select:

* **Ligand-Binding Domain:** A menu would display a list of known antibody fragments (e.g., scFv, VHH) that bind to the target ligand selected in Step 1\. The tool would display key "efficiency" metrics from the database, such as binding affinity (Kd​) and on-rate, to help the user make an informed choice.  
* **Transmembrane Domain:** The user could select from a library of transmembrane domains (e.g., from CD28, Notch) that are known to work well in synthetic receptors.  
* **Intracellular Domains:**  
  * **For the Protease Chain:** The user would choose a protease (e.g., TEV protease). The tool could show a metric for its cleavage efficiency.  
  * **For the Target Chain:** The user would select a cleavage site and a transcription factor (e.g., tTA, Cas9-based).

**c. Step 3: Define the Cellular Output** Finally, the user would specify the transcription factor that is released upon dimerization. This represents the "output" of the system. This TF could express a fluorescent reporter protein (like GFP) for research purposes or a therapeutic protein for a more advanced application.

### **2\. Backend & Logic**

The backend is where the magic happens. It would be responsible for storing the data and processing the user's selections.

**a. Database Management** The tool's most critical component is its database. This database would store all the information about the individual MESA components. It would be structured to hold data on:

* **Ligands:** Names and properties of molecules to be detected.  
* **Ligand-Binding Domains:** A comprehensive library of antibody fragments (scFv, VHH) and their specific binding partners. For each antibody, the database would store crucial "efficiency" metrics, such as:  
  * **Binding Affinity (**Kd​**)**: A measure of how tightly the antibody binds to the ligand. A lower Kd​ indicates stronger binding.  
  * **On-Rate (**kon​**)** and **Off-Rate (**koff​**)**: Metrics describing how fast the binding and unbinding occurs.  
  * **Specificity:** Information on potential off-target binding.  
  * https://github.com/bio-ontology-research-group/NanoDesigner  
* **Proteases and Cleavage Sites:** A list of known orthogonal protease-cleavage site pairs and their measured cleavage efficiencies.  
* **Transcription Factors:** A catalog of transcription factors and their cognate DNA binding sites.  
* **Transmembrane Domains:** Data on various transmembrane domains and their effects on receptor dimerization and stability.

**b. MESA Assembly Engine** When a user makes their selections, the backend would assemble the chosen components into a full MESA design. It would perform checks to ensure that the selected protease matches the cleavage site, and that the chosen transcription factor is compatible with the intended output gene.

### **3\. Output & Visualization**

After the user has completed their design, the tool would generate a detailed summary and a visual representation.

* **Interactive Diagram:** A dynamic SVG or canvas rendering of the MESA receptor, showing the two chains, the ligand binding, and the resulting cleavage and nuclear translocation.  
* **Performance Summary:** A detailed report that combines the efficiency metrics of the selected individual components to provide an estimated overall performance for the full MESA system. For example, it might provide a prediction for the signal-to-noise ratio or the sensitivity of the receptor.  
* **Design Recommendations:** The tool could suggest alternative components that might improve performance, such as a higher-affinity antibody or a more efficient protease. It could also flag potential issues with the chosen combination.

Citations:  
Dunbar, J., Krawczyk, K. et al (2014). Nucleic Acids Res. 42\. D1140-D1146  
