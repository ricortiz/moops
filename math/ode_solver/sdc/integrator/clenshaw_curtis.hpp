#ifndef CLENSHAW_CURTIS_HPP
#define CLENSHAW_CURTIS_HPP

// this file was generated by "mksmat.py"
#include<cassert>

namespace detail{

    template<int sdcnodes>
    inline long double clenshaw_curtis_quadrature_matrix(unsigned i, unsigned j);

    template<int sdcnodes>
    inline long double clenshaw_curtis_quadrature_node(unsigned i);

    template<int sdcnodes>
    inline long double clenshaw_curtis_quadrature_dt(unsigned i);

    template<>
    inline long double clenshaw_curtis_quadrature_matrix<2>(unsigned i, unsigned j)
    {
        static long double integration_matrix[1][2] = 
        {
            {0.5L,0.5L}
        };
        return integration_matrix[i][j];
    } // clenshaw_curtis_quadrature_matrix<2>(i,j)

    template<>
    inline long double clenshaw_curtis_quadrature_node<2>(unsigned i)
    {
        static long double integration_nodes[2] = {0.0L,1.0L};
        return integration_nodes[i];

    }// clenshaw_curtis_quadrature_nodes<2>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_dt<2>(unsigned i)
    {
        static long double integration_dt[1] = {1.0L};
        return integration_dt[i];

    }// clenshaw_curtis_quadrature_dt<2>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_matrix<3>(unsigned i, unsigned j)
    {
        static long double integration_matrix[2][3] = 
        {
            {0.20833333333333333333333333333333333L,0.33333333333333333333333333333333333L,-0.041666666666666666666666666666666667L},
            {-0.041666666666666666666666666666666667L,0.33333333333333333333333333333333333L,0.20833333333333333333333333333333333L}
        };
        return integration_matrix[i][j];
    } // clenshaw_curtis_quadrature_matrix<3>(i,j)

    template<>
    inline long double clenshaw_curtis_quadrature_node<3>(unsigned i)
    {
        static long double integration_nodes[3] = {0.0L,0.5L,1.0L};
        return integration_nodes[i];

    }// clenshaw_curtis_quadrature_nodes<3>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_dt<3>(unsigned i)
    {
        static long double integration_dt[2] = {0.5L,0.5L};
        return integration_dt[i];

    }// clenshaw_curtis_quadrature_dt<3>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_matrix<4>(unsigned i, unsigned j)
    {
        static long double integration_matrix[3][4] = 
        {
            {0.10243055555555555555555555555555556L,0.16319444444444444444444444444444444L,-0.024305555555555555555555555555555556L,0.0086805555555555555555555555555555557L},
            {-0.055555555555555555555555555555555556L,0.30555555555555555555555555555555556L,0.30555555555555555555555555555555556L,-0.055555555555555555555555555555555556L},
            {0.0086805555555555555555555555555555545L,-0.024305555555555555555555555555555554L,0.16319444444444444444444444444444444L,0.10243055555555555555555555555555556L}
        };
        return integration_matrix[i][j];
    } // clenshaw_curtis_quadrature_matrix<4>(i,j)

    template<>
    inline long double clenshaw_curtis_quadrature_node<4>(unsigned i)
    {
        static long double integration_nodes[4] = {0.0L,0.25L,0.75L,1.0L};
        return integration_nodes[i];

    }// clenshaw_curtis_quadrature_nodes<4>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_dt<4>(unsigned i)
    {
        static long double integration_dt[3] = {0.25L,0.5L,0.25L};
        return integration_dt[i];

    }// clenshaw_curtis_quadrature_dt<4>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_matrix<5>(unsigned i, unsigned j)
    {
        static long double integration_matrix[4][5] = 
        {
            {0.059701779686442458740014072701747484L,0.095031716019062009094954263719320677L,-0.012132034355964257320253308631454712L,0.0066433683707435685448487184562145472L,-0.002798220313557541259985927298252516L},
            {-0.043035113019775792073347406035080817L,0.21507831261090820533859016014022492L,0.21213203435596425732025330863145471L,-0.050086730334047116311726475649093474L,0.019464886980224207926652593964919183L},
            {0.019464886980224207926652593964919183L,-0.050086730334047116311726475649093474L,0.21213203435596425732025330863145471L,0.21507831261090820533859016014022492L,-0.043035113019775792073347406035080818L},
            {-0.0027982203135575412599859272982525158L,0.0066433683707435685448487184562145471L,-0.012132034355964257320253308631454713L,0.095031716019062009094954263719320676L,0.059701779686442458740014072701747484L}
        };
        return integration_matrix[i][j];
    } // clenshaw_curtis_quadrature_matrix<5>(i,j)

    template<>
    inline long double clenshaw_curtis_quadrature_node<5>(unsigned i)
    {
        static long double integration_nodes[5] = {0.0L,0.14644660940672623779957781894757548L,0.5L,0.85355339059327376220042218105242452L,1.0L};
        return integration_nodes[i];

    }// clenshaw_curtis_quadrature_nodes<5>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_dt<5>(unsigned i)
    {
        static long double integration_dt[4] = {0.14644660940672623779957781894757548L,0.35355339059327376220042218105242452L,0.35355339059327376220042218105242452L,0.14644660940672623779957781894757548L};
        return integration_dt[i];

    }// clenshaw_curtis_quadrature_dt<5>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_matrix<7>(unsigned i, unsigned j)
    {
        static long double integration_matrix[6][7] = 
        {
            {0.027217641773063450651486911537124147L,0.043301987114152472024853676078676587L,-0.0049824159104277846004735414275810216L,0.0023168848381701680174186430475824364L,-0.0015101936882055623782513192053587995L,0.0012035299857422600293949108337421896L,-0.00056013600471432712629086624065363049L},
            {-0.023646213201634879222915482965695576L,0.11130059290333706380765806959210959L,0.10900027305328492745761639857043816L,-0.02057085309213842198567261130155069L,0.011778050831062705235394176348215943L,-0.0089807131778349704650812596791315465L,0.0041315645761428985548622948120822023L},
            {0.017460317460317460317460317460317461L,-0.042997994092957230059923950470226231L,0.14915674603174603174603174603174603L,0.14841269841269841269841269841269841L,-0.034871031746031746031746031746031744L,0.023156724251687388790082680628956399L,-0.010317460317460317460317460317460319L},
            {-0.010317460317460317460317460317460317L,0.023156724251687388790082680628956401L,-0.034871031746031746031746031746031733L,0.14841269841269841269841269841269841L,0.14915674603174603174603174603174604L,-0.042997994092957230059923950470226249L,0.017460317460317460317460317460317457L},
            {0.004131564576142898554862294812082202L,-0.0089807131778349704650812596791315272L,0.011778050831062705235394176348215959L,-0.02057085309213842198567261130155069L,0.10900027305328492745761639857043817L,0.11130059290333706380765806959210962L,-0.023646213201634879222915482965695588L},
            {-0.00056013600471432712629086624065363037L,0.0012035299857422600293949108337421987L,-0.0015101936882055623782513192053587914L,0.0023168848381701680174186430475824367L,-0.0049824159104277846004735414275810177L,0.043301987114152472024853676078676598L,0.027217641773063450651486911537124144L}
        };
        return integration_matrix[i][j];
    } // clenshaw_curtis_quadrature_matrix<7>(i,j)

    template<>
    inline long double clenshaw_curtis_quadrature_node<7>(unsigned i)
    {
        static long double integration_nodes[7] = {0.0L,0.066987298107780676618138414623531908L,0.25L,0.5L,0.75L,0.93301270189221932338186158537646809L,1.0L};
        return integration_nodes[i];

    }// clenshaw_curtis_quadrature_nodes<7>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_dt<7>(unsigned i)
    {
        static long double integration_dt[6] = {0.066987298107780676618138414623531909L,0.18301270189221932338186158537646809L,0.25000000000000000000000000000000001L,0.25000000000000000000000000000000001L,0.18301270189221932338186158537646815L,0.066987298107780676618138414623531938L};
        return integration_dt[i];

    }// clenshaw_curtis_quadrature_dt<7>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_matrix<9>(unsigned i, unsigned j)
    {
        static long double integration_matrix[8][9] = 
        {
            {0.015446942589330503000508663449987255L,0.024571482764589171269473279957615346L,-0.0027292031693117551469555730885524011L,0.0012032717847525820515086196769462349L,-0.0007254858795644505687305583642067501L,0.0005193332763892971933218119544627156L,-0.00041852184693699231788565624639968958L,0.00037047163577776345415915451176587754L,-0.00017805741066949699949133655001274502L},
            {-0.014284667911499342383051696632943004L,0.065936015302512478001394667403264558L,0.064433133068355454270914110191314059L,-0.011352775377787173411387730111410457L,0.0059796718502256497349685233567425182L,-0.0040641423785444213585649548332606578L,0.0031968866471017310751071631737572791L,-0.002798077626495417682658972268551566L,0.0013403320885006576169483033670569956L},
            {0.012137064339490624190989302960854163L,-0.029504918979320072720126753003180839L,0.09681395624632992106397229434301373L,0.095937513323380720744747485105547128L,-0.020065348060354333061956683400867557L,0.011658919281530739014576110710998294L,-0.0086164639990830224068874253060961936L,0.0073388879192636753198885546661019445L,-0.0034879356605093758090106970391458451L},
            {-0.0093310850490678165544780158096444382L,0.020651758228406044235552779596616959L,-0.029751051205224818551352531429101463L,0.11352170202280769331051897659152622L,0.11322386050239154659413141682103022L,-0.026564892572284546796719478807455749L,0.016753803941309164552770158044604414L,-0.013456294636724564371715297182731607L,0.0062939149509321834455219841903555566L},
            {0.0062939149509321834455219841903555279L,-0.013456294636724564371715297182731863L,0.016753803941309164552770158044604307L,-0.02656489257228454679671947880745535L,0.11322386050239154659413141682102966L,0.11352170202280769331051897659152612L,-0.029751051205224818551352531429101765L,0.02065175822840604423555277959661683L,-0.0093310850490678165544780158096444724L},
            {-0.003487935660509375809010697039145936L,0.0073388879192636753198885546661022262L,-0.0086164639990830224068874253060961643L,0.011658919281530739014576110710999555L,-0.02006534806035433306195668340086967L,0.095937513323380720744747485105547664L,0.096813956246329921063972294343014024L,-0.029504918979320072720126753003179116L,0.012137064339490624190989302960854117L},
            {0.001340332088500657616948303367056873L,-0.0027980776264954176826589722685514887L,0.0031968866471017310751071631737572239L,-0.0040641423785444213585649548332590465L,0.0059796718502256497349685233567392369L,-0.01135277537778717341138773011140933L,0.064433133068355454270914110191313952L,0.065936015302512478001394667403264912L,-0.01428466791149934238305169663294307L},
            {-0.00017805741066949699949133655001280258L,0.00037047163577776345415915451176594706L,-0.00041852184693699231788565624639970263L,0.00051933327638929719332181195446349925L,-0.00072548587956445056873055836420848331L,0.0012032717847525820515086196769467959L,-0.0027292031693117551469555730885526981L,0.024571482764589171269473279957615764L,0.015446942589330503000508663449987241L}
        };
        return integration_matrix[i][j];
    } // clenshaw_curtis_quadrature_matrix<9>(i,j)

    template<>
    inline long double clenshaw_curtis_quadrature_node<9>(unsigned i)
    {
        static long double integration_nodes[9] = {0.0L,0.038060233744356621935908405301605857L,0.14644660940672623779957781894757548L,0.30865828381745511413577000798480057L,0.5L,0.69134171618254488586422999201519943L,0.85355339059327376220042218105242452L,0.96193976625564337806409159469839414L,1.0L};
        return integration_nodes[i];

    }// clenshaw_curtis_quadrature_nodes<9>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_dt<9>(unsigned i)
    {
        static long double integration_dt[8] = {0.038060233744356621935908405301605844L,0.10838637566236961586366941364596973L,0.16221167441072887633619218903722482L,0.19134171618254488586422999201520012L,0.191341716182544885864229992015199L,0.1622116744107288763361921890372267L,0.10838637566236961586366941364596926L,0.038060233744356621935908405301605561L};
        return integration_dt[i];

    }// clenshaw_curtis_quadrature_dt<9>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_matrix<13>(unsigned i, unsigned j)
    {
        static long double integration_matrix[12][13] = 
        {
            {0.0069091554849653485976588362077622526L,0.010989242416242542103495308538716176L,-0.0011905937159861112930031456228420269L,0.00050586775331751122553695563157662185L,-0.000289317317564458015095710482621544L,0.00019308756397099453995123507443877401L,-0.00014230004322922589287907593539263935L,0.00011266653440511920437336940473039781L,-0.000094356478176003952260917120595746046L,0.000082801230488540605442435753892059643L,-0.000075688805540651660506792727392403341L,0.000071811192051347009201509649963005952L,-0.000035288959479095846785608236682259685L},
            {-0.0066749341629101714066863242685191341L,0.030396588223050213267560446874821141L,0.029659853920918820447715258904875633L,-0.0049877287771198289944320968735405789L,0.0024868264863807902223483610623915058L,-0.001572263653641560264980786872320418L,0.0011276500473284478735131137613858219L,-0.0008791307108311391032262203309278112L,0.00072936564699233615951356770036551502L,-0.00063629124769892768322177498416313622L,0.0005795453466008628468070119538619972L,-0.00054878014828929640965866261618611594L,0.00026951028153427303775812017592582842L},
            {0.006222453314828841402579857769176388L,-0.014997336800523680380014089987559404L,0.047453516040898854179488908764151671L,0.046903603330951205816625313258253975L,-0.0091237430287923083899782645348084664L,0.0049050426361681303290490038120303824L,-0.0032840069800381201187887718150524584L,0.0024691673939039426973459085603950062L,-0.0020056874732367528344227089792523023L,0.0017273367551440028687935901237797713L,-0.001561116187178749786652367913877696L,0.001472073426435798439277611942103939L,-0.00072199112961560304186458667526943965L},
            {-0.0055825487627581444676782438342936816L,0.012248514613399648112978982202092485L,-0.017174612352449338531934650174149322L,0.061484467155043049483206356014774453L,0.061079112231604347811097242326666909L,-0.012876038180860002843026583928081246L,0.0072696119468938690931256889600149988L,-0.0050739083481009913960471375668420592L,0.0039610566760487922555416867711101437L,-0.0033336544537238069202043771781743213L,0.0029725064161469772087491874787810234L,-0.0027830120296569375821521706296811602L,0.0013618956816862999767662006101534232L},
            {0.0047988289444294263522028145956904454L,-0.010164328764809466064416915015790652L,0.012297827042327344081849698828302491L,-0.01854060596122353228638052635997527L,0.071401422420095010525502991872737772L,0.071150501600392918170809113619768289L,-0.015884552672618415066513871653698448L,0.0093233604845671461554252772409826957L,-0.0067449273081823324262206903781745044L,0.005452686601157022480408805587879419L,-0.0047543509238734100917846151955818931L,0.0044002314864929250869101278946753453L,-0.0021456155000150180922416298487578224L},
            {-0.0039247030703035522263286887215674014L,0.0081720467252998453280149368020613099L,-0.0092928445905953950257729987896073106L,0.011894065775572884586140270671764704L,-0.018811044048466638897131363501111248L,0.076496157746158138132379231468336947L,0.076411063199128941577040382180209423L,-0.017910796175610526321427815203115851L,0.010897805679810704054592318749804764L,-0.0081708557802157394895332592643908961L,0.0068572235799965688908157702647417861L,-0.0062283378636557416524551215681150469L,0.0030197413741408922181157557228842546L},
            {0.0030197413741408922181157557228690079L,-0.0062283378636557416524551215680131956L,0.0068572235799965688908157702647568534L,-0.0081708557802157394895332592644116153L,0.010897805679810704054592318749788754L,-0.017910796175610526321427815203079238L,0.076411063199128941577040382180159224L,0.07649615774615813813237923146835502L,-0.018811044048466638897131363501090182L,0.011894065775572884586140270671784853L,-0.0092928445905953950257729987895923569L,0.0081720467252998453280149368021801546L,-0.0039247030703035522263286887215785936L},
            {-0.0021456155000150180922416298487857921L,0.0044002314864929250869101278945722623L,-0.0047543509238734100917846151955362878L,0.0054526866011570224804088055877952197L,-0.0067449273081823324262206903782506287L,0.0093233604845671461554252772410928789L,-0.015884552672618415066513871653952636L,0.071150501600392918170809113619738106L,0.071401422420095010525502991872758661L,-0.018540605961223532286380526360021195L,0.012297827042327344081849698828323521L,-0.010164328764809466064416915016031138L,0.0047988289444294263522028145956780937L},
            {0.0013618956816862999767662006100784665L,-0.0027830120296569375821521706297467389L,0.0029725064161469772087491874789316373L,-0.0033336544537238069202043771782894153L,0.0039610566760487922555416867709118896L,-0.005073908348100991396047137566522229L,0.0072696119468938690931256889592972103L,-0.012876038180860002843026583928277023L,0.061079112231604347811097242326739681L,0.061484467155043049483206356014774579L,-0.017174612352449338531934650174042355L,0.012248514613399648112978982202091161L,-0.005582548762758144467678243834355527L},
            {-0.0007219911296156030418645866753842317L,0.0014720734264357984392776119418502971L,-0.0015611161871787497866523679135998132L,0.0017273367551440028687935901236579163L,-0.0020056874732367528344227089795975717L,0.0024691673939039426973459085609770186L,-0.003284006980038120118788771816363337L,0.0049050426361681303290490038116675984L,-0.0091237430287923083899782645346703992L,0.046903603330951205816625313258146296L,0.04745351604089885417948890876428229L,-0.014997336800523680380014089988176103L,0.0062224533148288414025798577691101818L},
            {0.00026951028153427303775812017580441063L,-0.0005487801482892964096586626164459996L,0.00057954534660086284680701195417899637L,-0.00063629124769892768322177498423826958L,0.0007293656469923361595135677000008142L,-0.00087913071083113910322622033027302751L,0.0011276500473284478735131137599323574L,-0.001572263653641560264980786872737065L,0.0024868264863807902223483610625209842L,-0.0049877287771198289944320968737508283L,0.029659853920918820447715258905045729L,0.030396588223050213267560446874919338L,-0.0066749341629101714066863242686069941L},
            {-0.000035288959479095846785608236735453971L,0.000071811192051347009201509649837470869L,-0.000075688805540651660506792727247961136L,0.000082801230488540605442435753868413282L,-0.000094356478176003952260917120756424796L,0.00011266653440511920437336940502842288L,-0.00014230004322922589287907593604902546L,0.00019308756397099453995123507424873741L,-0.00028931731756445801509571048256669423L,0.00050586775331751122553695563148166597L,-0.0011905937159861112930031456227510591L,0.010989242416242542103495308538716258L,0.0069091554849653485976588362077443305L}
        };
        return integration_matrix[i][j];
    } // clenshaw_curtis_quadrature_matrix<13>(i,j)

    template<>
    inline long double clenshaw_curtis_quadrature_node<13>(unsigned i)
    {
        static long double integration_nodes[13] = {0.0L,0.017037086855465856625128400135551316L,0.066987298107780676618138414623531908L,0.14644660940672623779957781894757548L,0.25L,0.37059047744873961882555058118797584L,0.5L,0.62940952255126038117444941881202416L,0.75L,0.85355339059327376220042218105242452L,0.93301270189221932338186158537646809L,0.98296291314453414337487159986444868L,1.0L};
        return integration_nodes[i];

    }// clenshaw_curtis_quadrature_nodes<13>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_dt<13>(unsigned i)
    {
        static long double integration_dt[12] = {0.017037086855465856625128400135552669L,0.049950211252314819993010014487970248L,0.079459311298945561181439404324071367L,0.10355339059327376220042218105237165L,0.12059047744873961882555058118805787L,0.12940952255126038117444941881189544L,0.12940952255126038117444941881212868L,0.12059047744873961882555058118738106L,0.10355339059327376220042218105159134L,0.079459311298945561181439404321900142L,0.049950211252314819993010014486350446L,0.01703708685546585662512840013481868L};
        return integration_dt[i];

    }// clenshaw_curtis_quadrature_dt<13>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_matrix<17>(unsigned i, unsigned j)
    {
        static long double integration_matrix[16][17] = 
        {
            {0.0038950713533365004403197865060139599L,0.0061950189195266359256232288113118945L,-0.00066539910616212945303957658549445881L,0.00027914146921208818526779299190424039L,-0.0001568058434915510076750746261045871L,0.00010222822224364891018823195945934367L,-0.000073173518185615681811119960470139492L,0.000055926780896448434865623559882090993L,-0.00004491704897099513340814448436960827L,0.00003752913014370816538674619077199704L,-0.000032404163988397377600953307716325469L,0.000028781857882482557987297372577130755L,-0.000026213453468915295781570529311795577L,0.000024424715021997985859391350846110261L,-0.000023246026549509653188479058193598994L,0.000022575157601877993595915235529322567L,-0.000011178646663499559680213494001473159L},
            {-0.0038204434440338403973082013975794384L,0.017316921142479962941960103133236032L,0.016888672857159681399769700772395515L,-0.0027952904841623698664531224236552115L,0.0013682761323604267196163439629736523L,-0.00084473816711935125258871087172573015L,0.00058820021053491252751988037201469617L,-0.00044245025376192048454185528973394778L,0.00035180550317430076074137348637080817L,-0.00029198988202139433473823937027355213L,0.00025096162999980954998348027144535519L,-0.00022218942381118846484735697691375668L,0.00020190078389873633244654729543232524L,-0.00018782561940001404305275619713340122L,0.00017857518795957715134703971692250204L,-0.00017331878325164164354650171569547248L,0.000085806555966159602691798602539024134L},
            {0.0036740555341312078364774231049163692L,-0.0088296820423001607743713851786221127L,0.02760131220634262594567217009112506L,0.027258499947660739397986306214245278L,-0.0051780567399383902087151644493913821L,0.0027157311565297458065982810197172734L,-0.0017636659398997718143190070994502465L,0.0012781301307240868330801932255018392L,-0.00099376334515644564078614390301505557L,0.00081293970700984328535223520862862087L,-0.00069189712309670764470479850195137527L,0.00060841524839560938121623359345263144L,-0.00055022785765450593121021271020149674L,0.0005101856719275726184986486665993719L,-0.00048401275080504845302040979838296085L,0.00046919076636915105046561329653243337L,-0.00023219446586879216352257689540419586L},
            {-0.003461533228910875256352581117375196L,0.0075734950696364471579975761360226035L,-0.010527383334684234095497316692387613L,0.036934547569809652301891609938216995L,0.036646177073402807432941605514369385L,-0.007477481695240087173982005894831289L,0.0040879650773713435642845002842954364L,-0.0027454382843440675377383418232174551L,0.0020473596946494784719425878010281269L,-0.0016325374419481017999141549286863021L,0.0013664113810838387587271074198661855L,-0.00118796965890721808025584648260342L,0.0010659672387379707744842728247398037L,-0.00098311353852704934162660821748315855L,0.0009294353351708490316609900808853724L,-0.00089920247039102261323880595965348195L,0.00044471677108912474364741888323871873L},
            {0.0031910436413785781854349756627981587L,-0.0067394515815607778906763915266486887L,0.0080830148842815524166369345593720479L,-0.012003180306336249323190332188158231L,0.044885014521793884513670929632164857L,0.044652322983471890869616128323876615L,-0.0095467345788427673235323842494207L,0.005376225547576385240411834657836714L,-0.0037041224178669088293560367532953558L,0.0028266524707197968127866644390531607L,-0.0023025738906227259903277983622395513L,0.001966661348444585298108844245125869L,-0.0017437387689107572135495962259288742L,0.0015953511376377657569643960148385776L,-0.0015004786322316872153435885539180062L,0.001447474083161506335917218724857847L,-0.00071520635862142181456502433818993844L},
            {-0.0028729815346908745063444027431102067L,0.005964130292526792786836723415653312L,-0.0067205170854548380375732416177849772L,0.0084693618327530608722888466720966844L,-0.013094351479827081547377143897050803L,0.051129027613631432176285793888126432L,0.050957527060864835124281018657654154L,-0.011285128662278881557302609065428492L,0.0065083762827865909730531582347100691L,-0.0045820920299614787784289531925275904L,0.0035686009718567383164589489452163603L,-0.002964782159524886775223213007100946L,0.0025820208584966077988910010118134228L,-0.0023347154987697593645324410428006656L,0.002179613158977884144892174380203992L,-0.0020939577594390406126758429154928932L,0.0010332684653091254936555972583099978L},
            {0.002519569857236488274094847894657716L,-0.005180069207987981939804018991550994L,0.0056480715986366009877741877808415269L,-0.0066204314790260149840688871570849317L,0.0086236143613959076150639805484186941L,-0.013725972310665783821852377376333673L,0.055419948686367632497633428849165461L,0.055314770531630100206426272104311863L,-0.012616738866415137479493526582993961L,0.007430950497204377982704014295068786L,-0.0053368420613876868238536240600804977L,0.004238503118391843280399753586540338L,-0.003591304443105587840869009291410657L,0.0031915561506760869942830598166602335L,-0.0029478062599608178084962020476339792L,0.0028154151442542365260508104837000772L,-0.001386680142763511725905152107301528L},
            {-0.0021443900215844394782826322237573846L,0.0043825140658723793225886598876000491L,-0.0046858812726632986623030896624416012L,0.0052796769362635046062918922177217844L,-0.0063665870228177279324400390096357818L,0.0085055665279018808962662306234228459L,-0.013860386311961289496991374521531859L,0.057589984111542563421019671498060994L,0.057554531653271748820526910715797239L,-0.013486540046668642454553854882530972L,0.0081049546881525303124951453812343499L,-0.0059377725207729523269573367906614263L,0.0048125488232447128482166582617945541L,-0.0041654122137955465745826257213395758L,0.003787195011558709507414902638057425L,-0.0035867013778955594062840519083935907L,0.0017618599784155605217173677783946838L},
            {0.0017618599784155605217173677771161854L,-0.0035867013778955594062840519306536587L,0.0037871950115587095074149026426961925L,-0.0041654122137955465745826257294268194L,0.0048125488232447128482166582519652036L,-0.0059377725207729523269573367838856264L,0.0081049546881525303124951453844130743L,-0.0134865400466686424545538548847681L,0.057554531653271748820526910714959227L,0.057589984111542563421019671497801585L,-0.01386038631196128949699137451613799L,0.0085055665279018808962662306218890865L,-0.006366587022817727932440039037387738L,0.0052796769362635046062918922067347827L,-0.0046858812726632986623030896565647319L,0.0043825140658723793225886598675419014L,-0.0021443900215844394782826322273078547L},
            {-0.0013866801427635117259051521015116181L,0.0028154151442542365260508104910614368L,-0.0029478062599608178084962020580166336L,0.0031915561506760869942830598360424488L,-0.0035913044431055878408690092379091604L,0.0042385031183918432803997536118724618L,-0.005336842061387686823853624053871102L,0.0074309504972043779827040142843920153L,-0.012616738866415137479493526584995385L,0.055314770531630100206426272101339473L,0.055419948686367632497633428850538083L,-0.013725972310665783821852377381201549L,0.0086236143613959076150639805683144813L,-0.0066204314790260149840688871363295099L,0.0056480715986366009877741877694528932L,-0.0051800692079879819398040189580520232L,0.0025195698572364882740948478962198312L},
            {0.0010332684653091254936555972660792549L,-0.0020939577594390406126758429270591335L,0.0021796131589778841448921743550500344L,-0.0023347154987697593645324409937306038L,0.0025820208584966077988910011460654956L,-0.0029647821595248867752232129335628606L,0.0035686009718567383164589489694143874L,-0.0045820920299614787784289532265036901L,0.0065083762827865909730531582373665921L,-0.011285128662278881557302609075020481L,0.050957527060864835124281018699703592L,0.051129027613631432176285793866713824L,-0.013094351479827081547377143988881948L,0.0084693618327530608722888466704210155L,-0.0067205170854548380375732416184792763L,0.0059641302925267927868367234036584558L,-0.0028729815346908745063444027538085085L},
            {-0.00071520635862142181456502432412164153L,0.0014474740831615063359172187778919007L,-0.001500478632231687215343588642605319L,0.0015953511376377657569643961855561037L,-0.0017437387689107572135495958403058398L,0.0019666613484445852981088444083017316L,-0.0023025738906227259903277983071247797L,0.002826652470719796812786664360326262L,-0.0037041224178669088293560367257408665L,0.0053762255475763852404118346347793957L,-0.0095467345788427673235323840965416803L,0.044652322983471890869616128252146704L,0.044885014521793884513670929691144145L,-0.012003180306336249323190332117146918L,0.0080830148842815524166369345273381073L,-0.0067394515815607778906763914514725652L,0.0031910436413785781854349756503647552L},
            {0.00044471677108912474364741888868330735L,-0.00089920247039102261323880580941538539L,0.00092943533517084903166098989535619806L,-0.00098311353852704934162660786092640315L,0.0010659672387379707744842735555181692L,-0.0011879696589072180802558462046954622L,0.001366411381083838758727107520463067L,-0.001632537441948101799914155061857991L,0.002047359694649478471942587885154337L,-0.0027454382843440675377383418670965733L,0.0040879650773713435642845005923362206L,-0.007477481695240087173982006056071016L,0.036646177073402807432941605439046892L,0.036934547569809652301891609932801182L,-0.010527383334684234095497316705948813L,0.0075734950696364471579975761619759844L,-0.0034615332289108752563525811521564993L},
            {-0.00023219446586879216352257690898223966L,0.00046919076636915105046561358092777085L,-0.00048401275080504845302041008100124172L,0.00051018567192757261849864920497409253L,-0.00055022785765450593121021167885992735L,0.00060841524839560938121623394728076994L,-0.00069189712309670764470479836827331848L,0.00081293970700984328535223504545411255L,-0.00099376334515644564078614375361407545L,0.0012781301307240868330801931630168261L,-0.0017636659398997718143190066484163196L,0.002715731156529745806598280773329298L,-0.0051780567399383902087151647697526221L,0.02725849994766073939798630625836582L,0.027601312206342625945672170026279152L,-0.0088296820423001607743713850550549162L,0.0036740555341312078364774230719544768L},
            {0.000085806555966159602691798574431774067L,-0.00017331878325164164354650141433409138L,0.00017857518795957715134703944073829231L,-0.00018782561940001404305275567474974496L,0.00020190078389873633244654825734044679L,-0.00022218942381118846484735666629637093L,0.00025096162999980954998348039301622766L,-0.00029198988202139433473823950677862951L,0.00035180550317430076074137364553789126L,-0.00044245025376192048454185534966131527L,0.00058820021053491252751988080447678792L,-0.00084473816711935125258871111500999679L,0.0013682761323604267196163436846562426L,-0.0027952904841623698664531223201772478L,0.016888672857159681399769700733808474L,0.017316921142479962941960103108484618L,-0.0038204434440338403973082014353410233L},
            {-0.000011178646663499559680213508865536299L,0.000022575157601877993595915369024657866L,-0.000023246026549509653188479175772513871L,0.000024424715021997985859391572239024674L,-0.000026213453468915295781570129157737542L,0.000028781857882482557987297497614635806L,-0.000032404163988397377600953257931587248L,0.000037529130143708165386746137597140786L,-0.000044917048970995133408144414211821921L,0.000055926780896448434865623534569103294L,-0.000073173518185615681811119778387437726L,0.00010222822224364891018823185573495191L,-0.00015680584349155100767507474500149323L,0.0002791414692120881852677930309943269L,-0.00066539910616212945303957658861861675L,0.0061950189195266359256232288045155508L,0.0038950713533365004403197864986105568L}
        };
        return integration_matrix[i][j];
    } // clenshaw_curtis_quadrature_matrix<17>(i,j)

    template<>
    inline long double clenshaw_curtis_quadrature_node<17>(unsigned i)
    {
        static long double integration_nodes[17] = {0.0L,0.0096073597983847754369088819328804816L,0.038060233744356621935908405301605857L,0.084265193848727381460605811191047122L,0.14644660940672623779957781894757548L,0.22221488349019888762858459302573356L,0.30865828381745511413577000798480057L,0.40245483899193586607585756576148888L,0.5L,0.59754516100806413392414243423851112L,0.69134171618254488586422999201519943L,0.77778511650980111237141540697426644L,0.85355339059327376220042218105242452L,0.91573480615127261853939418880895288L,0.96193976625564337806409159469839414L,0.99039264020161522456309111806711952L,1.0L};
        return integration_nodes[i];

    }// clenshaw_curtis_quadrature_nodes<17>(i)

    template<>
    inline long double clenshaw_curtis_quadrature_dt<17>(unsigned i)
    {
        static long double integration_dt[16] = {0.0096073597983847754369088819326341032L,0.028452873945971846498999523370619399L,0.046204960104370759524697405884300051L,0.062181415557998856338972007766424711L,0.075768274083472649829006774062124501L,0.08644340032725622650718541498248785L,0.093796555174480751940087557744974473L,0.097545161008064133924142434281791733L,0.097545161008064133924142434198984719L,0.093796555174480751940087557897346143L,0.086443400327256226507185415097426149L,0.075768274083472649829006774982789495L,0.062181415557998856338972009153167214L,0.046204960104370759524697407807627658L,0.028452873945971846498999525160142335L,0.0096073597983847754369088827029532042L};
        return integration_dt[i];

    }// clenshaw_curtis_quadrature_dt<17>(i)

} // namespace detail 

template<typename value_type, int sdc_nodes>
struct ClenshawCurtis
{
    inline value_type matrix(int i, int j)
    {
        assert(i < sdc_nodes-1 && j < sdc_nodes);
        return static_cast<value_type>(detail::clenshaw_curtis_quadrature_matrix<sdc_nodes>(i,j));
    }

    inline value_type node(int i)
    {
        assert(i < sdc_nodes);
        return static_cast<value_type>(detail::clenshaw_curtis_quadrature_node<sdc_nodes>(i));
    }

    inline value_type dt(int i)
    {
        assert(i < sdc_nodes-1);
        return static_cast<value_type>(detail::clenshaw_curtis_quadrature_dt<sdc_nodes>(i));
    }
};

#endif
