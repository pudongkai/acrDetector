# coding=utf-8


def classify(product):
    product = product.lower()
    product_type = {'conserved hypothetical': 0, 'hypothetical': 1, 'uncharacterised': 2,
                   'conserved putative': 0, 'putative': 1,
                   'helix-turn-helix': 3,
                   'phage': 4}
    types = ['conserved hypothetical', 'conserved putative', 'hypothetical', 'putative', 'uncharacterised',
             'helix-turn-helix', 'phage']
    for item in types:
        if item in product:
            return product_type[item]

    return 5