import unittest
# from scripts.entrez.fetch_gb_records import fetch_sources, fetch_gb_records
from fetch_gb_records import fetch_sources, fetch_gb_records

# Constants
EXPECTED_SINGLE_SOURCE = {'34577062': [{'authors': ['Sivendran,S. and Colman,R.F.'], 'title': 'Effect of a new non-cleavable substrate analog on wild-type and', 'journal': 'Protein Sci. 17 (7), 1162-1174 (2008)', 'year': 2008}, {'authors': ['Zhang,F.', 'Xu,Y.', 'Liu,P.', 'Fan,H.', 'Huang,X.', 'Sun,G.', 'Song,Y. and'], 'title': 'Association analyses of the interaction between the ADSS and ATM', 'journal': 'BMC Med. Genet. 9, 119 (2008)', 'year': 2008}, {'authors': ['Zhang,F.', 'Sham,P.C.', 'Fan,H.', 'Xu,Y.', 'Huang,X.', 'So,H.', 'Song,Y. and'], 'title': 'An association study of ADSS gene polymorphisms with schizophrenia', 'journal': 'Behav Brain Funct 4, 39 (2008)', 'year': 2008}, {'authors': ['Sun,H.', 'Li,N.', 'Wang,X.', 'Chen,T.', 'Shi,L.', 'Zhang,L.', 'Wang,J.', 'Wan,T.'], 'title': 'Molecular cloning and characterization of a novel muscle', 'journal': 'Mol. Cell. Biochem. 269 (1-2), 85-94 (2005)', 'year': 2005}, {'authors': ['Stepinski,J.', 'Pawelczyk,T.', 'Bizon,D. and Angielski,S.'], 'title': 'Purine nucleotide cycle in rat renal cortex and medulla under', 'journal': 'Kidney Int. 50 (4), 1195-1201 (1996)', 'year': 1996}, {
    'authors': ['Powell,S.M.', 'Zalkin,H. and Dixon,J.E.'], 'title': 'Cloning and characterization of the cDNA encoding human', 'journal': 'FEBS Lett. 303 (1), 4-10 (1992)', 'year': 1992}, {'authors': ['Zalkin,H. and Dixon,J.E.'], 'title': 'De novo purine nucleotide biosynthesis', 'journal': 'Prog. Nucleic Acid Res. Mol. Biol. 42, 259-287 (1992)', 'year': 1992}, {'authors': ['Lai,L.W.', 'Hart,I.M. and Patterson,D.'], 'title': 'A gene correcting the defect in the CHO mutant Ade -H, deficient in', 'journal': 'Genomics 9 (2), 322-328 (1991)', 'year': 1991}, {'authors': ['Ogawa,H.', 'Shiraki,H.', 'Matsuda,Y. and Nakagawa,H.'], 'title': 'Interaction of adenylosuccinate synthetase with F-actin', 'journal': 'Eur. J. Biochem. 85 (2), 331-337 (1978)', 'year': 1978}, {'authors': ['Van der Weyden,M.B. and Kelly,W.N.'], 'title': 'Human adenylosuccinate synthetase. Partial purification, kinetic', 'journal': 'J. Biol. Chem. 249 (22), 7282-7289 (1974)', 'year': 1974}]}

EXPECTED_MULTIPLE_SOURCES = {
    '34577062': [
        {
            'authors': ['Sivendran,S. and Colman,R.F.'],
            'title': 'Effect of a new non-cleavable substrate analog on wild-type and',
            'journal': 'Protein Sci. 17 (7), 1162-1174 (2008)',
            'year': 2008
        },
        {
            'authors': ['Zhang,F.', 'Xu,Y.', 'Liu,P.', 'Fan,H.', 'Huang,X.', 'Sun,G.', 'Song,Y. and'],
            'title': 'Association analyses of the interaction between the ADSS and ATM',
            'journal': 'BMC Med. Genet. 9, 119 (2008)',
            'year': 2008
        },
        {
            'authors': ['Zhang,F.', 'Sham,P.C.', 'Fan,H.', 'Xu,Y.', 'Huang,X.', 'So,H.', 'Song,Y. and'],
            'title': 'An association study of ADSS gene polymorphisms with schizophrenia',
            'journal': 'Behav Brain Funct 4, 39 (2008)',
            'year': 2008
        },
        {
            'authors': ['Sun,H.', 'Li,N.', 'Wang,X.', 'Chen,T.', 'Shi,L.', 'Zhang,L.', 'Wang,J.', 'Wan,T.'],
            'title': 'Molecular cloning and characterization of a novel muscle',
            'journal': 'Mol. Cell. Biochem. 269 (1-2), 85-94 (2005)',
            'year': 2005
        },
        {
            'authors': ['Stepinski,J.', 'Pawelczyk,T.', 'Bizon,D. and Angielski,S.'],
            'title': 'Purine nucleotide cycle in rat renal cortex and medulla under',
            'journal': 'Kidney Int. 50 (4), 1195-1201 (1996)',
            'year': 1996
        },
        {
            'authors': ['Powell,S.M.', 'Zalkin,H. and Dixon,J.E.'],
            'title': 'Cloning and characterization of the cDNA encoding human',
            'journal': 'FEBS Lett. 303 (1), 4-10 (1992)',
            'year': 1992
        },
        {'authors': ['Zalkin,H. and Dixon,J.E.'],
         'title': 'De novo purine nucleotide biosynthesis',
         'journal': 'Prog. Nucleic Acid Res. Mol. Biol. 42, 259-287 (1992)',
         'year': 1992
         },
        {
            'authors': ['Lai,L.W.', 'Hart,I.M. and Patterson,D.'],
            'title': 'A gene correcting the defect in the CHO mutant Ade -H, deficient in',
            'journal': 'Genomics 9 (2), 322-328 (1991)',
            'year': 1991
        },
        {
            'authors': ['Ogawa,H.', 'Shiraki,H.', 'Matsuda,Y. and Nakagawa,H.'],
            'title': 'Interaction of adenylosuccinate synthetase with F-actin',
            'journal': 'Eur. J. Biochem. 85 (2), 331-337 (1978)',
            'year': 1978
        },
        {
            'authors': ['Van der Weyden,M.B. and Kelly,W.N.'],
            'title': 'Human adenylosuccinate synthetase. Partial purification, kinetic',
            'journal': 'J. Biol. Chem. 249 (22), 7282-7289 (1974)',
            'year': 1974
        }
    ],
    '377581039': [{'authors': ['Arakaki,M.', 'Christin,P.A.', 'Nyffeler,R.', 'Lendel,A.', 'Eggli,U.,'], 'title': "Contemporaneous and recent radiations of the world's major", 'journal': 'Proc. Natl. Acad. Sci. U.S.A. 108 (20), 8379-8384 (2011)', 'year': 2011}, {'authors': ['Arakaki,M.', 'Christin,P.-A.', 'Nyffeler,R.', 'Eggli,U.', 'Ogburn,R.M.,'], 'title': 'Direct Submission', 'journal': 'Submitted (15-NOV-2010) Department of Ecology and Evolutionary'}]}

EXPECTED_RECORD_COUNT = '1686'
EXPECTED_RECORD_IDS = ['1066585321', '1393953329', '1519311736', '1519244926', '2462587281', 
                       '2462587279', '2462587277', '2462587276', '2462587274', '2194974489', 
                       '2194973393', '2194973003', '568815595', '568815587', '568815583', 
                       '2217341942', '2217341940', '2217341938', '2217341937', '2217341936', 
                       '2579155578', '251831106', '2234442025', '1111690406', '2032035867', 
                       '1969490283', '1969490269', '1969490255', '1969490241', '1969490227', 
                       '1969490213', '1969490199', '1969490185', '1969490171', '1969490157', 
                       '1969490143', '1969490129', '1969490115', '1969490101', '1969490087', 
                       '1969490073', '1969490059', '1969490045', '1969490031', '1969490017', 
                       '1969490003', '1969489989', '1969489975', '1969489961', '1969489947', 
                       '1969489933', '1969489919', '1002095861', '1875975052', '1875975038', 
                       '1875975024', '1875975010', '1875974996', '1875974982', '1875974968', 
                       '1875974954', '1875974940', '1875974926', '1875974912', '1875974898', 
                       '1875974884', '1875974870', '1875974856', '1875974842', '1875974828', 
                       '1875974814', '1875974800', '1875974786', '1875974772', '1875974758', 
                       '1875974744', '1875974730', '1875974716', '1875974702', '1875974688', 
                       '1875974674', '1875974660', '1875974646', '1875974632', '1875974618', 
                       '1875974604', '1875974590', '1875974576', '1875974562', '1875974548', 
                       '1875974534', '1875974520', '1875974506', '1875974492', '1875974478', 
                       '1875974464', '1875974450', '1875974436', '1875974422', '1875974408']

GI_NUM_1 = "34577062"
GI_NUM_2 = "377581039"
DATABASE = "nuccore"
SINGLE_ACCESSION = [GI_NUM_1]
MULTIPLE_ACCESSIONS = [GI_NUM_1, GI_NUM_2]
LOCUS = 'COI'
TAXID = "9606"


class TestFetchRecords(unittest.TestCase):

    # @patch('fetch_gb_records.fetch_metadata')
    def test_fetch_single_source(self):
        # pdb.set_trace()
        # mock_fetch_metadata.return_value = MOCK_METADATA

        result = fetch_sources(SINGLE_ACCESSION)
        print("ACTUAL SIGNLE SOURCES: ", result)
        print("EXPECTED SIGNLE SOURCE: ", EXPECTED_SINGLE_SOURCE)
        self.assertEqual(result, EXPECTED_SINGLE_SOURCE)

    def test_fetch_multiple_source(self):
        # pdb.set_trace()
        # mock_fetch_metadata.return_value = MOCK_METADATA
        result = fetch_sources(MULTIPLE_ACCESSIONS)
        print("ACTUAL SOURCES: ", result)
        print("EXPECTED SOURCES: ", EXPECTED_MULTIPLE_SOURCES)
        self.assertEqual(result, EXPECTED_MULTIPLE_SOURCES)

    # @patch('fetch_gb_records.fetch_')
    def test_fetch_gb_records_count(self):
        result = fetch_gb_records(LOCUS, TAXID, True)
        print("ACTUAL RECORD ID: ", result)
        print("EXPECTED RECORD ID: ", EXPECTED_RECORD_COUNT)
        self.assertEqual(result, EXPECTED_RECORD_COUNT)

    def test_fetch_gb_records_ids(self):
        result = fetch_gb_records(LOCUS, TAXID, False)
        print("ACTUAL RECORD ID: ", result)
        print("EXPECTED RECORD ID: ", EXPECTED_RECORD_IDS)
        self.assertEqual(result, EXPECTED_RECORD_IDS)
