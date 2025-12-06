class FrozenPipes(Exception):
    """
    Custom exception for WaterTAP developer errors.

    Raised when there's an issue that indicates a problem with the model
    implementation or configuration, rather than user input.

    Parameters
    ----------
    message : str
        Error message describing the issue
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"FrozenPipes: {self.message}"
