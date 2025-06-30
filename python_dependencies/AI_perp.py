import os
from openai import OpenAI

# ------------------------------------------------------------------
# Update these paths as needed to match your folder structure
# ------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))  # e.g., python_dependencies
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)               # one level above
API_KEY_FILE = os.path.join(PROJECT_ROOT, "api_keys", "api_key_perp.txt")
# ------------------------------------------------------------------

if not os.path.isfile(API_KEY_FILE):
    print(f"Error: {API_KEY_FILE} not found. Please create the file with your Perplexity API key.")
    exit(1)

with open(API_KEY_FILE, "r", encoding="utf-8") as f:
    api_key = f.read().strip()

# Instantiate the Perplexity API client using the OpenAI client-compatible interface
client = OpenAI(api_key=api_key, base_url="https://api.perplexity.ai")

# Define the maximum tokens to use (adjust as needed based on API limits)
MAX_TOKENS = 4096

print("------------------------------------------------------------------------------------------------------------------------------------")  
print("Hey there! I'm your AI assistant, here to help.")  
print("The future isn’t about man versus machine, but man with machine..!")  
print("Let's team up and create something great together!")  
print("------------------------------------------------------------------------------------------------------------------------------------")  

def query_perplexity(prompt, max_tokens=MAX_TOKENS, temperature=0.7, stream=False):
    """
    Query the Perplexity API with a given prompt.
    The model "sonar-pro" (top-tier model) is used.
    """
    messages = [
        {"role": "system", "content": "Be precise and concise. Provide ACS Nano-style citations with DOIs."},
        {"role": "user", "content": prompt}
    ]
    try:
        response = client.chat.completions.create(
            model="sonar-pro",
            messages=messages,
            max_tokens=max_tokens,
            temperature=temperature,
            stream=stream
        )
        # For non-streaming responses:
        if not stream:
            return response.choices[0].message.content.strip()
        else:
            # For streaming, accumulate and print the content as it arrives
            full_response = ""
            for part in response:
                chunk = part.choices[0].delta.content or ""
                full_response += chunk
                print(chunk, end="")
            return full_response
    except Exception as e:
        return f"API request failed: {e}"

def interactive_mode():
    """Option 1: Send a custom prompt to the Perplexity API."""
    print("\n--- Interactive Mode ---")
    user_input = input("Enter your prompt for Perplexity API: ")
    print("\nProcessing your request...\n")
    response = query_perplexity(user_input)  # now using MAX_TOKENS as default
    print("\nResponse from Perplexity API:\n")
    print(response)
    print("-" * 50)

def read_file_mode():
    """
    Option 2: List and display files from the MyAIAssistant folder.
    Update 'folder' below if you want to read from another location.
    """
    folder = "MyAIAssistant"
    if not os.path.isdir(folder):
        print(f"Directory '{folder}' does not exist.")
        return

    files = os.listdir(folder)
    if not files:
        print(f"No files found in '{folder}'.")
        return

    print("\nFiles available in MyAIAssistant:")
    for idx, fname in enumerate(files, 1):
        print(f"{idx}. {fname}")

    try:
        choice = int(input("Enter the number of the file you want to read: "))
        if 1 <= choice <= len(files):
            file_path = os.path.join(folder, files[choice - 1])
            with open(file_path, "r", encoding="utf-8") as file:
                content = file.read()
            print(f"\n--- Content of {files[choice - 1]} ---\n")
            print(content)
            print("-" * 50)
        else:
            print("Invalid file number.")
    except ValueError:
        print("Invalid input. Please enter a number.")

def quantum_espresso_analysis():
    """
    Option 3: List all regular files in the current folder (outputs)
    and let the user select one or more files (via comma-separated numbers)
    to be analyzed by the Perplexity API. The prompt instructs the AI to:
      - Analyze the Quantum Espresso output
      - Compare it with current literature
      - Provide a well-structured analysis
      - Include suggestions for improvement
      - Provide references in ACS Nano style with DOI
    """
    # List all regular files in the current folder
    all_files = [f for f in os.listdir(".") if os.path.isfile(f)]
    if not all_files:
        print("No files found in the current folder.")
        return

    print("\nFiles in the current folder:")
    for idx, fname in enumerate(all_files, 1):
        print(f"{idx}. {fname}")

    selection = input("Enter the numbers of the file(s) to analyze (comma-separated): ").strip()
    try:
        indices = [int(num.strip()) for num in selection.split(",")]
    except ValueError:
        print("Invalid input. Please enter numbers separated by commas.")
        return

    for idx in indices:
        if idx < 1 or idx > len(all_files):
            print(f"Invalid file number: {idx}. Skipping.")
            continue
        file_name = all_files[idx - 1]
        with open(file_name, "r", encoding="utf-8") as f:
            file_content = f.read()

        # Construct the prompt with detailed instructions for a structured scientific analysis
        prompt = (
            "You are an expert in computational materials science. "
            "Analyze the following Quantum Espresso output by comparing it with current scientific literature. "
            "Provide a well-structured response with the following sections:\n\n"
            "1. **Introduction**: Brief overview of the system and method.\n"
            "2. **Analysis**: Discuss the validity, implications, and any discrepancies with recent studies.\n"
            "3. **Detailed Suggestions**: Offer improvements to the simulation parameters, approaches, or methods, "
            "   referencing relevant literature.\n"
            "4. **References**: Provide references in ACS Nano style with DOI (e.g., 'ACS Nano, DOI: 10.1021/acsnano.xxxxxx').\n\n"
            f"Quantum Espresso Output from file '{file_name}':\n{file_content}\n\n"
            "Make sure your suggestions are backed by recent literature. Include relevant DOIs in the References section."
        )

        print(f"\n--- Analysis for file: {file_name} ---\n")
        analysis = query_perplexity(prompt, max_tokens=MAX_TOKENS)
        print(analysis)
        print("-" * 50)

def main():
    while True:
        print("\n=== Perplexity API Scientific Analysis Menu ===")
        print("1. Interactive mode (custom prompt)")
        print("2. Read a file from MyAIAssistant folder")
        print("3. Quantum Espresso analysis (analyze files in the current folder)")
        print("4. Exit")
        choice = input("Enter your choice (1-4): ").strip()

        if choice == "1":
            interactive_mode()
        elif choice == "2":
            read_file_mode()
        elif choice == "3":
            quantum_espresso_analysis()
        elif choice == "4":
            print("Exiting. Goodbye!")
            break
        else:
            print("Invalid choice. Please try again.")

if __name__ == "__main__":
    main()