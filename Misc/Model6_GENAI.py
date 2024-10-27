import os
from groq import Groq
import PyPDF2
import pandas as pd


api_key = "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"
client = Groq(api_key=api_key)

# Function to extract text from a fixed PDF file
def extract_text_from_pdf(pdf_path=r"C:\Users\Mishti mattu\OneDrive\Desktop\intel\pubchem_fingerprints.pdf"):
    with open(pdf_path, 'rb') as file:
        reader = PyPDF2.PdfReader(file)
        text = ""
        for page in reader.pages:
            text += page.extract_text()
    return text

# Function to extract text from a user-provided CSV file
def extract_text_from_csv(csv_path):
    df = pd.read_csv(csv_path)
    return df.to_string()

# Main function to handle the PDF, CSV, and custom prompt
def create_prompt_with_pdf_and_csv(csv_path, custom_prompt):
    # Extract text from the fixed PDF
    pdf_text = extract_text_from_pdf()

    # Extract text from the user-provided CSV
    csv_text = extract_text_from_csv(csv_path)

    # Combine the extracted PDF, CSV content, and the custom prompt
    combined_prompt = (
        f"PDF Content:\n{pdf_text}\n\n"
        f"CSV Content:\n{csv_text}\n\n"
        f"Custom Prompt:\n{custom_prompt}"
    )
    
    return combined_prompt

# Example usage of the above functions
def send_request_to_groq(csv_path, custom_prompt):
    # Read only the first 1000 rows of the CSV file
    df = pd.read_csv(csv_path, nrows=100)
    
    # Convert the DataFrame to a string
    csv_content = df.to_csv(index=False)
    
    # Your existing code to prepare the messages
    messages = [
        {"role": "system", "content": custom_prompt},
        {"role": "user", "content": f"Here's the CSV content:\n\n{csv_content}"}
    ]
    
    # Your existing code to send the request
    chat_completion = client.chat.completions.create(
        model="llama3-8b-8192",
        messages=messages,
        temperature=0.5,
        max_tokens=1024,
        top_p=1,
        stream=False,
        stop=None
    )
    
    return chat_completion

# User provides the CSV file and a custom prompt
csv_path = r"acetylcholinesterase_04_bioactivity_data_3class_pIC50.csv"  # Replace with the actual CSV path
custom_prompt = """
 I need to generate a new molecule that incorporates the following fingerprints from the provided bioactivity and fingerprint data files.
     528,3,493,193,601,335,390,813,308,180,614,750,758,559,259
     The pdf file contains fingerprint information 
     Please generate a new molecule that incorporates all the 15 fingerprints listed above.
     Return the canonical SMILES in of the newly generated molecule.
     The generated molecule should be optimized for biological activity based on the data provided.
     Ensure that all 15 fingerprints are present in the generated molecule's structure and the new molecule should be undiscovered
     ONLY RETURN THE SMILE OF THE NEWLY GENERATED MOLECULE"""
     #Evaluate whether this molecule can realistically exist based on its chemical and physical properties.
    
      

# Send the request with the CSV file and custom prompt
response = send_request_to_groq(csv_path, custom_prompt)  # Store the response
print("Custom Prompt:", custom_prompt)  # Print the custom prompt
print("Response:", response)  # Print the response


'''
# Define a list of multiple prompts
prompts = [
    "Can you give me the latest news?",
]

# Iterate through the list of prompts and get responses
for prompt in prompts:
    chat_completion = client.chat.completions.create(
        messages=[
            {
                "role": "user",
                "content": prompt
            }
        ],
        model="llama3-8b-8192",
    )

    # Print the response for each prompt
    print(f"Prompt: {prompt}")
    print(f"Response: {chat_completion.choices[0].message.content}\n")
'''

