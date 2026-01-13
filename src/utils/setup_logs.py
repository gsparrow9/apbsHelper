import logging
import warnings
def setup_log(level=logging.INFO):
    if not logging.getLogger().handlers:
        logging.basicConfig(
            level=level,
            format='[%(asctime)s] (%(levelname)s) %(name)s: %(message)s'
        )

        logging.getLogger("MDAnalysis").setLevel(logging.WARNING)
        warnings.filterwarnings(
            'ignore',
            module='MDAnalysis',
            category=UserWarning
        )

