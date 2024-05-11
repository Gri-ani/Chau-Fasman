# Chau-Fasman
### Protein Secondary Structure Prediction

Proteins are essential molecules in biological systems, performing a multitude of functions vital for life. Understanding their structure is crucial for unraveling their biological roles. Proteins often fold into intricate shapes, comprising different structural elements such as alpha helices and beta strands. Predicting these structural elements from a protein sequence can provide valuable insights into its function.

In this project, we develop a Python program to predict the secondary structure of a protein sequence. We utilize propensity values, which represent the likelihood of certain amino acids to form specific structural motifs, such as alpha helices or beta strands. The program employs a sliding window approach to analyze segments of the protein sequence, identifying regions with high propensity values indicative of alpha helices or beta strands.

#### Function Definitions

The code defines several functions to facilitate the analysis of protein sequences:

- **Alpha and Beta Properties**: Initializes dictionaries (`alpha_props` and `beta_props`) storing propensity values for amino acids to form alpha helices and beta strands, respectively.
  
- **Finding Alpha Helices**: Identifies potential alpha helices by examining sets of six amino acids at a time, considering regions with high propensity values.
  
- **Finding Beta Strands**: Identifies potential beta strands using a similar sliding window approach but with sets of five amino acids.
  
- **Extending Helices and Strands**: Extends identified helices or strands by adding adjacent amino acids, stopping when the propensity values fall below a threshold.
  
- **Printing Functions**: Functions to print out various aspects of the protein sequence and its secondary structures in a readable format.

#### Main Program

The program processes a provided protein sequence, predicts alpha helices and beta strands, and prints them out. It also identifies conflicting regions where both alpha helices and beta strands might occur simultaneously.

This Python program offers a valuable tool for researchers and bioinformaticians to quickly analyze protein sequences and gain insights into their structural characteristics. By leveraging predefined propensity values, it provides a convenient way to predict secondary structure elements, aiding in the understanding of protein function and behavior.
