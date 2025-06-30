import curses
import time

# Game-like text lines
TEXT_LINES = [
    "=============================================================",
    "WELCOME to the Neuromorphic Materials Calculator 2025!",
    "=============================================================",
    "Developed by Dr. Ing. Santiago D. Barrionuevo",
    "Under the supervision of Dra. Myriam H. Aguirre",
    "-------------------------------------------------------------",
    "This project is funded by the EU project MELON as part of",
    "HORIZON-2020: Marie Curie Research and Innovation Staff Exchange Action",
    "Empowering students and researchers with first-principles tools",
    "for simulating materials in neuromorphic applications (and more..!).",
    "-------------------------------------------------------------",
    "Memristive and multiferroic materials for logic units in nanoelectronics",
    "",
    "MELON MISSION",
    "The most simplistic computational model of a neuron is an 'on-off' switch,",
    "with a '0' representing a resting state and a '1' representing an axon firing",
    "an action potential. While this lends itself well to conventional digital",
    "electronics and silicon-based transistors, it does not represent the incredible",
    "natural 'state' space of a real neuron. When it comes to realising the potential",
    "of a brain-like processing system, novel materials are needed.",
    "The EU-funded MELON project has created an expert consortium of academic",
    "institutions and an SME to explore novel materials with history-dependent",
    "conductivity to emulate neuronal connectivity. Together with materials capable",
    "of multivalued logic and interconnects, the team plans to deliver the building",
    "blocks of tomorrow's emergent computing circuits.",
    "",
    "This program is crucial in simulating the materials that compose the heart",
    "of these neuromorphic devices, helping researchers explore their electronic",
    "and spintronic properties through first-principles calculations.",
    "-------------------------------------------------------------",
    "",
    "How it works:",
    "- Simulation & AI Assistance for Material Science:",
    "    Leverage our tool to simulate material properties and assist your research projects.",
    "",
    "- Literature Review:",
    "    Our AI searches relevant literature to support your project with up-to-date", 
    "    information and proper citations without hallucinations.",
    "",
    "- Multiple Modes:",
    "    (1) Existing Materials:",
    "        If you already know the material you want to simulate", 
    "        simply enter its Materials Project ID.",
    "    (2) Exploring Options:",
    "        If you're uncertain about the ID or are exploring potential materials",
    "        use our AI engine to find the ideal candidate.",
    "    (3) Pushing the Frontier of Science:",
    "        For experts creating new materials or pushing research boundaries",
    "        our AI assists in designing your ideal candidate.",
    "",
    "- Seamless Workflow:",
    "    The program automatically handles the following:",
    "        - AI search engine and assitance to the user.",
    "        - Generation of Quantum ESPRESSO input files.",
    "        - Execution of SCF, NSCF, bands, DOS, and EELS calculations.",
    "        - Processing and visualization of the resulting data.",
    "",
    "- Future Enhancements:",
    "    Stay tuned—more functionalities are coming soon!",
    "============================================================="
]

# Function to display text with animation
def draw_text(stdscr):
    stdscr.clear()
    stdscr.nodelay(False)  # Allow user input (blocking mode)
    stdscr.keypad(True)
    
    # Colors
    curses.start_color()
    curses.init_pair(1, curses.COLOR_CYAN, curses.COLOR_BLACK)    # Title
    curses.init_pair(2, curses.COLOR_GREEN, curses.COLOR_BLACK)   # Body
    curses.init_pair(3, curses.COLOR_YELLOW, curses.COLOR_BLACK)  # Highlighted section
    curses.init_pair(4, curses.COLOR_RED, curses.COLOR_BLACK)     # Red text
    
    height, width = stdscr.getmaxyx()
    
    for i, line in enumerate(TEXT_LINES):
        color = curses.color_pair(2)  # Default: Green

        if "WELCOME" in line or "MELON MISSION" in line:
            color = curses.color_pair(1)  # Cyan for titles
        elif "How it works:" in line or line.startswith("===") or line.startswith("---"):
            color = curses.color_pair(3)  # Yellow for separators and section headers
        elif any(sub in line for sub in ["- Literature Review:", "- Multiple Modes:", "- Seamless Workflow:", "- Future Enhancements:"]):
            color = curses.color_pair(4)  # Red for specific sections
        
        # Center the text
        x = max(0, (width // 2) - (len(line) // 2))
        stdscr.addstr(i + 2, x, line, color)
        stdscr.refresh()
        time.sleep(0.08)  # Simulate typing effect
    
    # Wait for user to press ENTER or SPACEBAR
    prompt = "[ Press ENTER or SPACE to continue ]"
    stdscr.addstr(height - 2, max(0, (width // 2) - (len(prompt) // 2)), prompt, curses.A_BOLD)
    stdscr.refresh()
    
    while True:
        key = stdscr.getch()
        if key in [curses.KEY_ENTER, 10, 13, 32]:  # 10, 13 = ENTER, 32 = SPACE
            break  # Exit loop when ENTER or SPACE is pressed

# Main function to start curses UI
def main():
    curses.wrapper(draw_text)

if __name__ == "__main__":
    main()
