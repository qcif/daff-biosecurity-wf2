"""Generic utility functions."""


def deduplicate(sequence, key=None):
    """Remove duplicates from a sequence while preserving order.

    If key() func is provided, uniqueness will be determined by key(element).
    If no key is provided, the element itself will be used to determine
    uniqueness.
    """
    if key:
        seen = set()
        return [
            item for item in sequence
            if key(item)
            and key(item) not in seen
            and not seen.add(key(item))
        ]
    return list(dict.fromkeys(sequence))
