import os
import timeit
import predict_interactions.interaction_clear_3did_kbdock as ic3k
import predict_interactions.process_tables as pt
import predict_interactions.filtering as filtering
import predict_interactions.pvalue as pv

result_address = 'resultdata/'
source_address = 'sourcedata/'


def main(params: dict):
    threshold_ignore = params.get('ignore_threshold', False)
    redo_similarity = params.get('redo_similarities', True)
    coefficents = params.get('coefficients', None)
    iterations = params.get('iterations', 10_000)
    sources = [x for x in os.listdir(source_address) if x.startswith('source')]
    print("Sources:", sources)

    start = timeit.default_timer()
    # This part gets all the domain-protein information (what proteins are associated with which domains etc.)
    seqDom, seqpdbchain, pdbchainDom = pt.read_chain_dom(source_address)
    # # pickle.dump((seqDom, seqpdbchain, pdbchainDom), open('../pickles/seqDom_seqpdbchain_pdbchainDom.pickle', 'wb'))
    # # seqDom, seqpdbchain, pdbchainDom = pickle.load(open('../pickles/seqDom_seqpdbchain_pdbchainDom.pickle', 'rb'))
    # # print("Loading files from pickle took:", round(timeit.default_timer() - start, 1), "seconds")
    #
    # This calculates all the similarity scores for each source
    for i in sources:
        pt.similarity_calculator_interaction(i, 'pfam', seqDom, seqpdbchain, pdbchainDom,
                                             source_address, result_address, redo=redo_similarity)

    # # sifts_reader_process('sifts', 'pfam')
    # This is a cleanup of the domain interaction sources
    ic3k.clean_3did_kbdock_domine_downloaded_files(source_address, result_address)

    # This function creates random wrong associations (retaining node degree) for the ddi inference
    filtering.create_wrong_assocations(sources, source_address, result_address)

    # This function assigns the interactions. It also does all the "hyperparameter optimization"
    filtering.assign_interaction(sources, result_address,
                                 continue_flag=threshold_ignore, coefficients=coefficents, iteration_option=iterations)

    ic3k.kbdock_union_3did(source_address, result_address)

    # This part calculates and sorts all the p-values for every domain-domain interaction
    for i in sources:
        pv.pvalue_calculation(i, seqDom, pdbchainDom, source_address, result_address)
    pv.accumulate_pvalues(sources, result_address)
    pv.gold_silver_bronze(result_address)

    # pv.one_to_one()
    print("Took:", timeit.default_timer() - start)
    # checkOneOnedomain()
    # verify_go_for_each_pair()
