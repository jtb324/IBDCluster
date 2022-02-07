
class UnSupportedIBDProgram(Exception):
    """Error that will be raised if the user doesn't provide a supported ibd program"""
    def __init__(self, program: str, message: str) -> None:
        self.program: str = program
        self.message: str = message
        super().__init__(message)

def check_ibd_program(program: str) -> str:
    """Callback function that will check if the ibd program provided is hapibd or ilash
    
    Parameters
    
    program : str
        string that has the value either 'hapibd' or 'ilash'

    Returns

    str
        returns the program if an error is not raise
    """
    if program.lower() not in ["hapibd", "ilash"]:
        raise UnSupportedIBDProgram(
            program, 
            f"The provided ibd program {program} is not supported. Supported values are 'hapibd' and 'ilash'."
            )

    return program.lower()