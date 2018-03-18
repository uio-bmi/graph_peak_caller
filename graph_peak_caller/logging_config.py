import logging
import sys
def set_logging_config():
    class InfoFilter(logging.Filter):
        def filter(self, rec):
            return rec.levelno in (logging.WARNING, logging.DEBUG, logging.INFO)

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    h1 = logging.StreamHandler(sys.stdout)
    h1.setLevel(logging.DEBUG)
    h1.addFilter(InfoFilter())
    formatter = logging.Formatter('%(asctime)s, %(levelname)s: %(message)s')
    h1.setFormatter(formatter)
    h2 = logging.StreamHandler()
    h2.setLevel(logging.WARNING)
    h2.setFormatter(formatter)

    logger.addHandler(h1)
    logger.addHandler(h2)
    logger.propagate = False

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s, %(levelname)s: %(message)s",
        handlers=[h1, h2]
    )