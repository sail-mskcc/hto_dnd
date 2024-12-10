import logging

# setup logging to the console in the format:
# <timestamp> <log level> - <message>
def get_logger(logger="hto_dnd", level=1):
    # Set logger level
    log = logging.getLogger(logger)
    level = {
        0: logging.WARNING,
        1: logging.INFO,
        2: logging.DEBUG,
    }.get(level, logging.INFO)
    log.setLevel(level)

    # Clear existing handlers
    if log.hasHandlers():
        log.handlers.clear()

    formatter = logging.Formatter("%(asctime)s HTO-DND | %(levelname)s: %(message)s")
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    log.addHandler(ch)

    return log