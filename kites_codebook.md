# kites.csv

One row per paragraph returned by a search of eHRAF World Cultures for "kite OR kiting" in November 2026 (582 rows, including paragraphs not about kites), and normally one row per society per survey for the six comparative surveys (250 rows). A 0 in any of the use columns means the record does not mention that use, not that it was absent.

pid: record identifier. For eHRAF, the page and row of the search results (p1_001 is page 1, row 1). For the surveys, tyl_ (Tylor 1880), bal_ (Balfour 1913), pli_ (Plischke 1922), lau_ (Laufer 1928), cha_ (Chadwick 1931) or ane_ (Anell 1955) and a number.
source: eHRAF, or the survey the record comes from.
society: the society the record describes.
locality: for the surveys, the place as the survey names it.
owc, culture: eHRAF Outline of World Cultures code and culture name.
region, subregion, subsistence: eHRAF classification of the culture.
ea_id, sccs_id: Ethnographic Atlas and Standard Cross-Cultural Sample society, with D-PLACE identifiers.
glottocode: Glottolog 5.3 language.
latitude, longitude: location of the society.
author, year, title, section, doc_page: the source and page. For eHRAF these identify the paragraph in eHRAF World Cultures.
ocms: eHRAF Outline of Cultural Materials subjects of the paragraph.
text: for eHRAF, the paragraph as returned, except that the 46 paragraphs over 500 words, among them whole indexes, glossaries and word lists, are cut to the passages around each mention of kites; for the surveys, the key quotation (in German for Plischke). In both, [...] marks words left out.
duplicate: 1 if the same paragraph already appears for the same culture (none do).
relevant: 1 if the record is about kites, the flying device on a line, that the society makes, flies, names, sells, believes in, tells of or forbids; 0 otherwise.
sense: what "kite" refers to: kite where relevant = 1, otherwise bird (the bird of prey), vernacular (a word in another language spelled kite), name, shape (something kite-shaped) or other.
narrative: 1 if the information comes from a myth, tale, song, proverb or other narrative rather than description.
play: flown as a pastime or game.
contest: kite fighting, or contests of height, duration or beauty.
gambling: wagers on kite contests.
fishing: kites used to catch fish.
ritual: kites in religious or magical practice.
divination: the flight or fall of a kite read as an omen.
festival: flown at a named festival or in a particular season.
myth: kites in myth, legend, tale or song.
sound: kites fitted with hummers, whistles, organs or rattles.
utility: practical uses other than fishing, such as signalling or warfare.
prohibition: rules restricting kite flying.
commerce: kites made, bought or sold for money.
construction: how kites are made, their materials, size or shape.
origin: a claim about where or when the society's kites came from.
denial: kites, or a kind of kite, said to be absent or no longer flown.
other: about kites but none of the above.
fliers: who flies the kites (children, boys, girls, youths, men, women, adults, elders, chiefs, priests or everyone), separated by semicolons.
materials: materials of the kite, its line or its lure, separated by semicolons.
form: shape of the kite (bird, human, animal, fish, insect, geometric, box or other).
term: vernacular word for kite given in the record.
denied: what the record says is absent, where denial = 1.
note: what the record says about kites in a line, or for relevant = 0, why it is not about kites.
scope: regional where a survey's statement about a whole island group or region is coded for a society within it.
link_note: the reason for any change to the society's links.
lineage: Glottolog 5.3 classification of the glottocode, from family to language, separated by slashes.

Records are coded for the people they describe. A people described at length in another people's file is coded for that people (Blackwood's account of Petats in the Tinputz file), while passing comparisons with other peoples, and tales or reports drawn from other sources about them, are not coded.

ritual2 (made in kites.R, not in the file): ritual with cha_15 (East Ambae), cha_16 (Malays) and cha_18 (Balinese) set to 0, since these rest only on a kite-song or a plate caption reported by Chadwick (1931). Used for the strict ritual sensitivity.
