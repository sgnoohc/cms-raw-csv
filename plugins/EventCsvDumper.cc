// The unpacked readout of CMS events as CSV files, one per subsystem: each row is one channel that was read out, with
// its address, its raw values (every time sample) and the global position of the channel in cm. Plus the raw bytes
// of every FED. For each file, extent/<same name> gives the outline of every channel, row by row: the four corners of
// its sensitive area. Also inventory/: every readout element of the detector, read out or not. Runs after RawToDigi;
// writes every event it sees, each into its own directory <dir>/run<run>_evt<event>/.
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/ESDataFrame.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

namespace {
  // the outline of a channel: the four corners of its sensitive area, in order around its edge
  using Box = std::optional<std::array<GlobalPoint, 4>>;
  template <class D>
  Box quad(const D& det, const LocalPoint& a, const LocalPoint& b, const LocalPoint& c, const LocalPoint& d) {
    return std::array<GlobalPoint, 4>{{det.toGlobal(a), det.toGlobal(b), det.toGlobal(c), det.toGlobal(d)}};
  }
  // strip positions s0 to s1 (in strip units from the edge of the strip plane) over the whole strip length
  template <class D, class T>
  Box strips(const D& det, const T& tp, float s0, float s1) {
    return quad(det, tp.localPosition(MeasurementPoint(s0, -0.5)), tp.localPosition(MeasurementPoint(s1, -0.5)),
                tp.localPosition(MeasurementPoint(s1, 0.5)), tp.localPosition(MeasurementPoint(s0, 0.5)));
  }
  // a calorimeter cell's face toward the interaction point: of its two faces (corners 0-3 and 4-7), the one nearer
  // the origin, with its corners put in order around its centre
  Box face(const CaloCellGeometry* cell) {
    if (!cell) return std::nullopt;
    const auto& cv = cell->getCorners();
    double d[2] = {0, 0};
    for (int i = 0; i < 8; ++i) d[i / 4] += cv[i].mag();
    const int o = d[1] < d[0] ? 4 : 0;
    std::array<GlobalPoint, 4> f{{cv[o], cv[o + 1], cv[o + 2], cv[o + 3]}};
    const GlobalPoint c((f[0].x() + f[1].x() + f[2].x() + f[3].x()) / 4, (f[0].y() + f[1].y() + f[2].y() + f[3].y()) / 4,
                        (f[0].z() + f[1].z() + f[2].z() + f[3].z()) / 4);
    const GlobalVector e1 = (f[0] - c).unit(), e2 = GlobalVector(c.x(), c.y(), c.z()).unit().cross(e1);
    std::sort(f.begin(), f.end(), [&](const GlobalPoint& a, const GlobalPoint& b) {
      return std::atan2((a - c).dot(e2), (a - c).dot(e1)) < std::atan2((b - c).dot(e2), (b - c).dot(e1));
    });
    return f;
  }
  // a CSV file that opens with '#' lines describing it and its columns (pandas: read_csv(path, comment='#')), and its
  // companion extent/<name> with the outline of each row's channel. xyz() ends a row: the global point in cm (or three
  // empty fields when there is none) and the outline (or twelve empty fields)
  struct Csv {
    std::ofstream o, e;
    Csv(const std::string& dir, const std::string& name, const std::string& head, const std::vector<std::string>& notes,
        const std::string& outline)
        : o(dir + "/" + name), e(dir + "/extent/" + name) {
      for (const auto& n : notes) o << "# " << n << '\n';
      o << std::fixed << std::setprecision(2) << head << '\n';
      e << "# " << notes[0] << '\n'
        << "# The outline of each channel in ../" << name << ": one row per row of that file and in the same order.\n"
        << "# " << outline << '\n'
        << "# x0_cm ... z3_cm: its four corners in global CMS coordinates in cm (x toward the LHC centre; y up; z along the beam) in order around its edge\n";
      e << std::fixed << std::setprecision(4) << "x0_cm,y0_cm,z0_cm,x1_cm,y1_cm,z1_cm,x2_cm,y2_cm,z2_cm,x3_cm,y3_cm,z3_cm\n";
    }
    template <class T>
    Csv& operator<<(const T& x) { o << x << ','; return *this; }
    void xyz(const GlobalPoint* g, const Box& b) {
      if (g) o << g->x() << ',' << g->y() << ',' << g->z() << '\n';
      else o << ",,\n";
      if (b)
        for (int i = 0; i < 4; ++i) e << (*b)[i].x() << ',' << (*b)[i].y() << ',' << (*b)[i].z() << (i < 3 ? ',' : '\n');
      else e << ",,,,,,,,,,,\n";
    }
  };
  std::string cols(const char* name, int n) {
    std::string s;
    for (int i = 0; i < n; ++i) s += std::string(name) + std::to_string(i) + ',';
    return s;
  }
  // a plain CSV with '#' lines on top, for the channel inventory
  struct Inv {
    std::ofstream o;
    Inv(const std::string& path, const std::vector<std::string>& notes, const std::string& head) : o(path) {
      for (const auto& n : notes) o << "# " << n << '\n';
      o << std::fixed << std::setprecision(3) << head << '\n';
    }
    template <class T>
    Inv& operator<<(const T& x) { o << x << ','; return *this; }
    Inv& p(const GlobalPoint& g) { o << g.x() << ',' << g.y() << ',' << g.z() << ','; return *this; }
    Inv& v(const GlobalVector& g) { o << g.x() << ',' << g.y() << ',' << g.z() << ','; return *this; }
    void end() { o.seekp(-1, std::ios_base::cur); o << '\n'; }
  };
  const char* fedSystem(int id) {
    const std::vector<std::tuple<const char*, int, int>> R = {
        {"pixel", 1200, 1349}, {"strip", 50, 489}, {"ecal_preshower", 520, 575}, {"ecal", 600, 670}, {"hcal", 700, 731},
        {"hcal", 1100, 1199}, {"csc", 750, 757}, {"csc", 830, 869}, {"dt", 770, 779}, {"dt", 1369, 1371}, {"rpc", 790, 795},
        {"gem", 1467, 1478}, {"ctpps", 577, 587}, {"ctpps", 1462, 1466}, {"tcds", 1024, 1099}, {"l1_trigger", 1350, 1368},
        {"l1_trigger", 1372, 1409}, {"daq_metadata", 1022, 1022}};
    for (const auto& [n, a, b] : R)
      if (a <= id && id <= b) return n;
    return "other";
  }
}  // namespace

class EventCsvDumper : public edm::one::EDAnalyzer<> {
public:
  explicit EventCsvDumper(const edm::ParameterSet& c)
      : dir_(c.getParameter<std::string>("dir")),
        rawT_(consumes(edm::InputTag("rawDataCollector"))),
        pixT_(consumes(edm::InputTag("siPixelDigis"))),
        stripT_(consumes(edm::InputTag("siStripDigis", "ZeroSuppressed"))),
        ebT_(consumes(edm::InputTag("ecalDigis", "ebDigis"))),
        eeT_(consumes(edm::InputTag("ecalDigis", "eeDigis"))),
        esT_(consumes(edm::InputTag("ecalPreshowerDigis"))),
        hbheT_(consumes(edm::InputTag("hcalDigis"))),
        hfT_(consumes(edm::InputTag("hcalDigis"))),
        hoT_(consumes(edm::InputTag("hcalDigis"))),
        dtT_(consumes(edm::InputTag("muonDTDigis"))),
        cscsT_(consumes(edm::InputTag("muonCSCDigis", "MuonCSCStripDigi"))),
        cscwT_(consumes(edm::InputTag("muonCSCDigis", "MuonCSCWireDigi"))),
        csccT_(consumes(edm::InputTag("muonCSCDigis", "MuonCSCComparatorDigi"))),
        rpcT_(consumes(edm::InputTag("muonRPCDigis"))),
        gemT_(consumes(edm::InputTag("muonGEMDigis"))),
        tkT_(esConsumes()), topoT_(esConsumes()), caloT_(esConsumes()), dtGT_(esConsumes()), cscGT_(esConsumes()),
        rpcGT_(esConsumes()), gemGT_(esConsumes()) {}

  void analyze(const edm::Event& e, const edm::EventSetup& es) override {
    const std::string dir = dir_ + "/run" + std::to_string(e.id().run()) + "_evt" + std::to_string(e.id().event());
    std::filesystem::create_directories(dir + "/extent");
    const auto& tk = es.getData(tkT_);
    const auto& topo = es.getData(topoT_);
    const auto& cg = es.getData(caloT_);
    const auto f = [&dir](const char* name) { return dir + "/" + name; };
    // the first comment line of every file, and the last one, which explains x_cm,y_cm,z_cm
    const std::string evt = "CMS run " + std::to_string(e.id().run()) + " lumi section " +
                            std::to_string(e.luminosityBlock()) + " event " + std::to_string(e.id().event()) +
                            ". Raw readout from the RAW data tier unpacked by CMSSW RawToDigi; values are raw codes before any calibration.";
    const auto pos = [](const std::string& what) {
      return "x_cm y_cm z_cm: " + what + " in global CMS coordinates in cm (x toward the LHC centre; y up; z along the beam)";
    };

    {  // tracker: one row per pixel or strip above threshold
      Csv p(dir, "pixel.csv", "part,layer,detid,row,col,adc,x_cm,y_cm,z_cm",
            {evt, "Silicon pixel tracker. One row per pixel read out (above the zero-suppression threshold).",
             "part: BPix = barrel pixel; FPix = forward pixel disks",
             "layer: BPix layer 1-4 (r = 2.9/6.8/10.9/16.0 cm) or FPix disk 1-3 (the sign of z gives the side)",
             "detid: CMSSW DetId of the sensor module (32-bit integer)",
             "row: pixel row in the module (0-159; 100 um pitch)",
             "col: pixel column in the module (0-415; 150 um pitch)",
             "adc: collected charge as an 8-bit code (0-255). One number per pixel; no time samples",
             pos("centre of the pixel")},
            "Pixel: the pixel cell (100 x 150 um; the pixels at the edges of a readout chip are larger).");
      for (const auto& ds : e.get(pixT_)) {
        const DetId id(ds.detId());
        const auto* du = dynamic_cast<const PixelGeomDetUnit*>(tk.idToDetUnit(id));
        for (const auto& d : ds) {
          p << (id.subdetId() == 1 ? "BPix" : "FPix") << topo.layer(id) << id.rawId() << d.row() << d.column() << d.adc();
          const GlobalPoint g = du ? du->surface().toGlobal(du->specificTopology().localPosition(MeasurementPoint(d.row() + 0.5, d.column() + 0.5))) : GlobalPoint();
          Box ol;
          if (du) {
            const auto at = [&](float r, float c) { return du->specificTopology().localPosition(MeasurementPoint(r, c)); };
            ol = quad(*du, at(d.row(), d.column()), at(d.row() + 1, d.column()), at(d.row() + 1, d.column() + 1), at(d.row(), d.column() + 1));
          }
          p.xyz(du ? &g : nullptr, ol);
        }
      }
      Csv s(dir, "strip.csv", "part,layer,stereo,detid,strip,adc,x_cm,y_cm,z_cm",
            {evt, "Silicon strip tracker. One row per strip read out (zero-suppressed by the readout boards).",
             "part: TIB/TOB = inner/outer barrel; TID = inner disks; TEC = endcaps",
             "layer: barrel layer (TIB 1-4; TOB 1-6) or disk (TID 1-3; TEC 1-9)",
             "stereo: 1 = the stereo sensor of a double-sided module (rotated by 100 mrad); 0 otherwise",
             "detid: CMSSW DetId of the sensor (32-bit integer)",
             "strip: strip number in the sensor (0-511 or 0-767)",
             "adc: charge code. 0-253 linear; 254 = 254 to 1022; 255 = saturated. One number per strip; no time samples",
             pos("middle of the strip")},
            "Strip: the strip over its whole length and one pitch wide (a strip measures one coordinate only).");
      const char* part[] = {"", "", "", "TIB", "TID", "TOB", "TEC"};
      for (const auto& ds : e.get(stripT_)) {
        const DetId id(ds.detId());
        const auto* du = dynamic_cast<const StripGeomDetUnit*>(tk.idToDetUnit(id));
        for (const auto& d : ds) {
          s << part[id.subdetId()] << topo.layer(id) << int(topo.isStereo(id)) << id.rawId() << d.strip() << d.adc();
          const GlobalPoint g = du ? du->surface().toGlobal(du->specificTopology().localPosition(d.strip() + 0.5)) : GlobalPoint();
          s.xyz(du ? &g : nullptr, du ? strips(*du, du->specificTopology(), d.strip(), d.strip() + 1) : Box());
        }
      }
    }
    // a calorimeter cell's reference point (into g; nullptr when it has no geometry) and its outline (into ol)
    const auto calo = [&cg](const DetId& id, GlobalPoint& g, Box& ol) -> const GlobalPoint* {
      const auto cell = cg.getGeometry(id);
      if (cell) g = cell->getPosition();
      ol = face(cell.get());
      return cell ? &g : nullptr;
    };
    const std::string facing = "the face of the cell toward the interaction point.";
    {  // ECAL: 10 samples of (12-bit ADC, gain); the gain is written as the factor 12, 6 or 1
      const int gain[] = {0, 12, 6, 1};
      const std::string sr = "One row per crystal read out. Selective readout keeps whole regions around energy deposits so many rows are pedestal only.";
      const std::string adc = "adc0-adc9: 10 time samples 25 ns apart; 12-bit ADC (0-4095). The pedestal is about 200 at gain 12; an in-time pulse peaks near sample 5";
      const std::string gn = "gain0-gain9: gain of each sample: 12 (normal) or 6 or 1 (large signals). The amplitude is roughly (adc - pedestal) x 12/gain; the pedestal differs per gain";
      Csv b(dir, "ecal_barrel.csv", "ieta,iphi," + cols("adc", 10) + cols("gain", 10) + "x_cm,y_cm,z_cm",
            {evt, "ECAL barrel (lead tungstate crystals). " + sr, "ieta: crystal index in eta (-85 to -1 and 1 to 85)",
             "iphi: crystal index in phi (1-360)", adc, gn, pos("centre of the crystal's front face")},
            "Crystal: " + facing);
      const auto& eb = e.get(ebT_);
      for (size_t i = 0; i < eb.size(); ++i) {
        const EcalDataFrame df(eb[i]);
        const EBDetId id(df.id());
        b << id.ieta() << id.iphi();
        for (int k = 0; k < df.size(); ++k) b << df.sample(k).adc();
        for (int k = 0; k < df.size(); ++k) b << gain[df.sample(k).gainId()];
        GlobalPoint g;
        Box ol;
        b.xyz(calo(id, g, ol), ol);
      }
      Csv n(dir, "ecal_endcap.csv", "zside,ix,iy," + cols("adc", 10) + cols("gain", 10) + "x_cm,y_cm,z_cm",
            {evt, "ECAL endcaps (lead tungstate crystals). " + sr, "zside: -1/+1 = endcap at negative/positive z",
             "ix iy: crystal position on the endcap's x-y grid (1-100 each)", adc, gn, pos("centre of the crystal's front face")},
            "Crystal: " + facing);
      const auto& ee = e.get(eeT_);
      for (size_t i = 0; i < ee.size(); ++i) {
        const EcalDataFrame df(ee[i]);
        const EEDetId id(df.id());
        n << id.zside() << id.ix() << id.iy();
        for (int k = 0; k < df.size(); ++k) n << df.sample(k).adc();
        for (int k = 0; k < df.size(); ++k) n << gain[df.sample(k).gainId()];
        GlobalPoint g;
        Box ol;
        n.xyz(calo(id, g, ol), ol);
      }
      Csv p(dir, "ecal_preshower.csv", "zside,plane,six,siy,strip," + cols("adc", 3) + "x_cm,y_cm,z_cm",
            {evt, "ECAL preshower (silicon strip sensors in front of the ECAL endcaps). One row per strip read out (zero-suppressed).",
             "zside: -1/+1 = negative/positive z",
             "plane: 1 = front plane; 2 = rear plane (its strips are perpendicular to those of plane 1)",
             "six siy: sensor position on the x-y grid (1-40 each)",
             "strip: strip number in the sensor (1-32)",
             "adc0-adc2: 3 time samples 25 ns apart. Pedestal-subtracted so values can be negative",
             pos("centre of the strip")},
            "Strip: the 1.9 mm x 6.1 cm silicon strip; " + facing);
      const auto& esd = e.get(esT_);
      for (size_t i = 0; i < esd.size(); ++i) {
        const ESDataFrame df(esd[i]);
        const ESDetId id(df.id());
        p << id.zside() << id.plane() << id.six() << id.siy() << id.strip();
        for (int k = 0; k < df.size(); ++k) p << df.sample(k).adc();
        GlobalPoint g;
        Box ol;
        p.xyz(calo(id, g, ol), ol);
      }
    }
    {  // HCAL: QIE codes per time sample; soi = index of the sample of the triggered crossing
      const std::string soi = "soi: index of the sample from the triggered bunch crossing (sample of interest; counted from 0)";
      Csv b(dir, "hcal_hbhe.csv", "part,ieta,iphi,depth,soi," + cols("adc", 8) + cols("tdc", 8) + "x_cm,y_cm,z_cm",
            {evt, "HCAL barrel (HB) and endcap (HE): brass and plastic scintillator read by SiPMs. One row per channel read out (zero-suppressed). Calibration and monitoring channels are left out.",
             "part: HB = barrel; HE = endcap",
             "ieta: tower index in eta (HB |ieta| 1-16; HE |ieta| 16-29)",
             "iphi: tower index in phi (1-72)",
             "depth: longitudinal segment counted from the inside (HB 1-4; HE 1-7)",
             soi,
             "adc0-adc7: 8 time samples 25 ns apart. 8-bit nonlinear QIE11 code; converting it to charge needs the lookup table from the conditions database",
             "tdc0-tdc7: arrival time in each sample. HE: 6-bit in 0.5 ns steps (0-49 valid; 62/63 = no valid time). HB: 2-bit code (3 = no valid time)",
             pos("front face of the cell")},
            "Cell: " + facing);
      const auto& hbhe = e.get(hbheT_);
      for (auto it = hbhe.begin(); it != hbhe.end(); ++it) {
        const QIE11DataFrame q(*it);
        const HcalDetId id(q.detid());
        if (id.subdet() != HcalBarrel && id.subdet() != HcalEndcap) continue;  // calibration and monitoring channels
        b << (id.subdet() == HcalBarrel ? "HB" : "HE") << id.ieta() << id.iphi() << id.depth() << q.presamples();
        for (int k = 0; k < q.samples(); ++k) b << q[k].adc();
        for (int k = 0; k < q.samples(); ++k) b << q[k].tdc();
        GlobalPoint g;
        Box ol;
        b.xyz(calo(id, g, ol), ol);
      }
      Csv h(dir, "hcal_hf.csv", "ieta,iphi,depth,soi," + cols("adc", 3) + cols("tdc_rise", 3) + cols("tdc_fall", 3) + "x_cm,y_cm,z_cm",
            {evt, "HCAL forward (HF): steel with quartz fibres read by PMTs. Every channel is read out in every event (no zero suppression). Calibration channels are left out.",
             "ieta: tower index in eta (|ieta| 29-41)",
             "iphi: tower index in phi (odd numbers 1-71; every 4th from 3 for |ieta| 40-41)",
             "depth: 1 and 3 = long fibres; 2 and 4 = short fibres. Each PMT has two anodes: depth 3 is the second anode of depth 1 and depth 4 of depth 2",
             soi,
             "adc0-adc2: 3 time samples 25 ns apart. 8-bit nonlinear QIE10 code",
             "tdc_rise0-tdc_rise2: rising-edge time in each sample; 6-bit in 0.5 ns steps (0-49 valid; 62/63 = no valid time)",
             "tdc_fall0-tdc_fall2: falling-edge time code in each sample (5-bit)",
             pos("front face of the cell")},
            "Cell: " + facing);
      const auto& hf = e.get(hfT_);
      for (auto it = hf.begin(); it != hf.end(); ++it) {
        const QIE10DataFrame q(*it);
        const HcalDetId id(q.detid());
        if (id.subdet() != HcalForward) continue;
        h << id.ieta() << id.iphi() << id.depth() << q.presamples();
        for (int k = 0; k < q.samples(); ++k) h << q[k].adc();
        for (int k = 0; k < q.samples(); ++k) h << q[k].le_tdc();
        for (int k = 0; k < q.samples(); ++k) h << q[k].te_tdc();
        GlobalPoint g;
        Box ol;
        h.xyz(calo(id, g, ol), ol);
      }
      Csv o(dir, "hcal_ho.csv", "ieta,iphi,depth,soi," + cols("adc", 10) + "x_cm,y_cm,z_cm",
            {evt, "HCAL outer (HO): scintillator outside the solenoid read by SiPMs. One row per channel read out (zero-suppressed).",
             "ieta: tower index in eta (-15 to -1 and 1 to 15)",
             "iphi: tower index in phi (1-72)",
             "depth: always 4 (HO is a single layer)",
             soi,
             "adc0-adc9: 10 time samples 25 ns apart; 7-bit nonlinear QIE8 code (0-127)",
             pos("front face of the cell")},
            "Tile: " + facing);
      for (const auto& df : e.get(hoT_)) {
        o << df.id().ieta() << df.id().iphi() << df.id().depth() << df.presamples();
        for (int k = 0; k < df.size(); ++k) o << df.sample(k).adc();
        GlobalPoint g;
        Box ol;
        o.xyz(calo(df.id(), g, ol), ol);
      }
    }
    {  // muon system
      const auto& dtg = es.getData(dtGT_);
      Csv d(dir, "dt.csv", "wheel,station,sector,superlayer,layer,wire,hit,tdc_counts,time_ns,x_cm,y_cm,z_cm",
            {evt, "Drift tubes (DT; muon barrel). One row per wire hit. Only a time is recorded; there is no amplitude.",
             "wheel: -2 to +2 (wheel 0 is at z = 0)",
             "station: 1-4 = MB1-MB4 from the inside out",
             "sector: 1-12 in phi (30 deg each); 13 and 14 are the second MB4 chambers of sectors 4 and 10",
             "superlayer: 1 and 3 measure r-phi; 2 measures z (MB4 has no superlayer 2)",
             "layer: 1-4 within the superlayer",
             "wire: wire number in the layer",
             "hit: index of the hit when a wire fired more than once (0 1 2 ...)",
             "tdc_counts: raw TDC time (1 count = 25/32 ns)",
             "time_ns: tdc_counts x 25/32. It still includes the trigger latency and per-chamber offsets so it is not yet the drift time",
             pos("midpoint of the wire")},
            "Drift cell: 4.2 cm wide around the wire and as long as the wire (the hit is somewhere in the cell; its time says how far from the wire).");
      const auto& dt = e.get(dtT_);
      for (auto r = dt.begin(); r != dt.end(); ++r) {
        const DTLayerId& id = (*r).first;
        const DTLayer* l = dtg.layer(id);
        for (auto x = (*r).second.first; x != (*r).second.second; ++x) {
          d << id.wheel() << id.station() << id.sector() << id.superlayer() << id.layer() << x->wire() << x->number() << x->countsTDC() << x->time();
          const GlobalPoint g = l ? l->toGlobal(LocalPoint(l->specificTopology().wirePosition(x->wire()), 0, 0)) : GlobalPoint();
          Box ol;
          if (l) {
            const auto& tp = l->specificTopology();
            const float w = tp.wirePosition(x->wire()), hw = tp.cellWidth() / 2, hl = tp.cellLenght() / 2;
            ol = quad(*l, LocalPoint(w - hw, -hl, 0), LocalPoint(w + hw, -hl, 0), LocalPoint(w + hw, hl, 0), LocalPoint(w - hw, hl, 0));
          }
          d.xyz(l ? &g : nullptr, ol);
        }
      }
      const auto& cscg = es.getData(cscGT_);
      const auto bins = [](const std::vector<int>& v) {
        std::string s;
        for (int b : v) s += (s.empty() ? "" : " ") + std::to_string(b);
        return s;
      };
      const std::string csc = "Cathode strip chambers (CSC; muon endcaps): ";
      const std::vector<std::string> cscId = {
          "zside: +1/-1 = endcap at positive/negative z", "station: 1-4 = ME1-ME4",
          "ring: 1-3 in station 1 and 1-2 in stations 2-4 (ring 1 is innermost)",
          "chamber: chamber number in phi (1-36; 1-18 in ring 1 of stations 2-4)", "layer: 1-6"};
      const auto cscNotes = [&](const std::string& what, std::vector<std::string> rest) {
        std::vector<std::string> v{evt, csc + what};
        v.insert(v.end(), cscId.begin(), cscId.end());
        v.insert(v.end(), rest.begin(), rest.end());
        return v;
      };
      Csv s(dir, "csc_strips.csv", "zside,station,ring,chamber,layer,strip," + cols("adc", 8) + "x_cm,y_cm,z_cm",
            cscNotes("cathode strips. Whole front-end boards are read out around a hit so many rows are pedestal only. "
                     "ME1/1 is split here: ring 1 = ME1/1b and ring 4 = ME1/1a (its high-eta end).",
                     {"strip: strip number in the layer (1-80; ME1/1b 1-64; ME1/1a 1-48)",
                      "adc0-adc7: 8 time samples 50 ns apart; 12-bit ADC. The pedestal is about 500-1000 and differs per strip",
                      pos("centre of the strip")}),
            "Strip: the cathode strip over its whole length and one strip wide.");
      const auto& cs = e.get(cscsT_);
      for (auto r = cs.begin(); r != cs.end(); ++r) {
        const CSCDetId& id = (*r).first;
        const CSCLayer* l = cscg.layer(id);
        for (auto x = (*r).second.first; x != (*r).second.second; ++x) {
          s << (id.endcap() == 1 ? 1 : -1) << id.station() << id.ring() << id.chamber() << id.layer() << x->getStrip();
          for (int a : x->getADCCounts()) s << a;
          // the middle of the strip; CSCLayer::centerOfStrip is on the layer centre line, off the strips of ME1/1a and ME1/1b
          const GlobalPoint g = l ? l->toGlobal(l->geometry()->topology()->localPosition(MeasurementPoint(x->getStrip() - 0.5, 0))) : GlobalPoint();
          s.xyz(l ? &g : nullptr, l ? strips(*l, *l->geometry()->topology(), x->getStrip() - 1, x->getStrip()) : Box());
        }
      }
      Csv w(dir, "csc_wires.csv", "zside,station,ring,chamber,layer,wiregroup,timebins_on,x_cm,y_cm,z_cm",
            cscNotes("anode wire groups. One row per wire group that fired.",
                     {"wiregroup: wire group number in the layer",
                      "timebins_on: the 25 ns time bins (0-15) in which the wire group fired; space-separated",
                      pos("centre of the wire group")}),
            "Wire group: the band of wires read out together, as long as its middle wire.");
      const auto& cw = e.get(cscwT_);
      for (auto r = cw.begin(); r != cw.end(); ++r) {
        const CSCDetId& id = (*r).first;
        const CSCLayer* l = cscg.layer(id);
        for (auto x = (*r).second.first; x != (*r).second.second; ++x) {
          w << (id.endcap() == 1 ? 1 : -1) << id.station() << id.ring() << id.chamber() << id.layer() << x->getWireGroup() << bins(x->getTimeBinsOn());
          const GlobalPoint g = l ? l->centerOfWireGroup(x->getWireGroup()) : GlobalPoint();
          Box ol;
          if (l) {
            const auto* lg = l->geometry();
            const int wg = x->getWireGroup();
            const LocalPoint m = lg->localCenterOfWireGroup(wg);
            const float a = lg->wireAngle(), hl = lg->lengthOfWireGroup(wg) / 2, hw = lg->numberOfWiresPerGroup(wg) * lg->wirePitch() / 2;
            const LocalVector u(std::cos(a) * hl, std::sin(a) * hl, 0), v(-std::sin(a) * hw, std::cos(a) * hw, 0);
            ol = quad(*l, m - u - v, m + u - v, m + u + v, m - u + v);
          }
          w.xyz(l ? &g : nullptr, ol);
        }
      }
      Csv c(dir, "csc_comparators.csv", "zside,station,ring,chamber,layer,strip,half,timebins_on,x_cm,y_cm,z_cm",
            cscNotes("comparator hits. The strip comparators find which half-strip a muon crossed; they feed the trigger.",
                     {"strip: strip number in the layer. ME1/1 (ring 1) counts 1-112: 1-64 are ME1/1b strips and 65-112 are ME1/1a strips 1-48",
                      "half: 0 = left half of the strip; 1 = right half",
                      "timebins_on: the 25 ns time bins (0-15) in which the comparator fired; space-separated",
                      pos("centre of the strip")}),
            "Half-strip: the half of the strip the comparator picked, over the whole strip length.");
      const auto& cc = e.get(csccT_);
      for (auto r = cc.begin(); r != cc.end(); ++r) {
        const CSCDetId& id = (*r).first;
        for (auto x = (*r).second.first; x != (*r).second.second; ++x) {
          c << (id.endcap() == 1 ? 1 : -1) << id.station() << id.ring() << id.chamber() << id.layer() << x->getStrip() << x->getComparator() << bins(x->getTimeBinsOn());
          // in data, ME1/1a comparators come in ring 1 as strips 65-112: they are strips 1-48 of the ring-4 (ME1/1a) layer
          const bool me1a = id.station() == 1 && id.ring() == 1 && x->getStrip() > 64;
          const CSCLayer* l = cscg.layer(me1a ? CSCDetId(id.endcap(), 1, 4, id.chamber(), id.layer()) : id);
          const int strip = me1a ? x->getStrip() - 64 : x->getStrip();
          const GlobalPoint g = l ? l->toGlobal(l->geometry()->topology()->localPosition(MeasurementPoint(strip - 0.5, 0))) : GlobalPoint();
          const float s0 = strip - 1 + 0.5f * x->getComparator();
          c.xyz(l ? &g : nullptr, l ? strips(*l, *l->geometry()->topology(), s0, s0 + 0.5f) : Box());
        }
      }
      const auto& rpcg = es.getData(rpcGT_);
      const std::string bx = "bx: bunch crossing of the hit relative to the triggered one (0)";
      Csv p(dir, "rpc.csv", "region,ring,station,sector,layer,subsector,roll,strip,bx,x_cm,y_cm,z_cm",
            {evt, "Resistive plate chambers (RPC; muon barrel and endcaps). One row per strip hit. Hit or no hit; no amplitude.",
             "region: 0 = barrel; -1/+1 = endcap at negative/positive z",
             "ring: the wheel in the barrel (-2 to +2); the ring in the endcap (2-3)",
             "station: 1-4",
             "sector: 1-12 in the barrel (30 deg each); 1-6 in the endcap (60 deg each)",
             "layer: 1 = inner and 2 = outer chamber in barrel stations 1-2; 1 elsewhere",
             "subsector: chamber within the sector (1-6 in the endcap at 10 deg each)",
             "roll: eta partition of the chamber (1-3)",
             "strip: strip number in the roll",
             bx,
             pos("centre of the strip")},
            "Strip: the strip over its whole length and one strip wide.");
      const auto& rpc = e.get(rpcT_);
      for (auto r = rpc.begin(); r != rpc.end(); ++r) {
        const RPCDetId& id = (*r).first;
        const RPCRoll* roll = rpcg.roll(id);
        for (auto x = (*r).second.first; x != (*r).second.second; ++x) {
          p << id.region() << id.ring() << id.station() << id.sector() << id.layer() << id.subsector() << id.roll() << x->strip() << x->bx();
          const GlobalPoint g = roll ? roll->toGlobal(roll->centreOfStrip(x->strip())) : GlobalPoint();
          p.xyz(roll ? &g : nullptr, roll ? strips(*roll, roll->specificTopology(), x->strip() - 1, x->strip()) : Box());
        }
      }
      const auto& gemg = es.getData(gemGT_);
      Csv m(dir, "gem.csv", "region,station,ring,chamber,layer,roll,strip,bx,x_cm,y_cm,z_cm",
            {evt, "Gas electron multipliers (GEM; muon endcaps: GE1/1 and GE2/1). One row per strip hit. Hit or no hit; no amplitude.",
             "region: -1/+1 = endcap at negative/positive z",
             "station: 1 = GE1/1; 2 = GE2/1",
             "ring: always 1",
             "chamber: chamber number in phi (1-36 for GE1/1 at 10 deg each; 1-18 for GE2/1 at 20 deg each)",
             "layer: 1-2 (each station is a pair of chambers)",
             "roll: eta partition of the chamber (1-8 in GE1/1)",
             "strip: strip number in the partition (0-383)",
             bx,
             pos("centre of the strip")},
            "Strip: the strip over its whole length and one strip wide.");
      const auto& gem = e.get(gemT_);
      for (auto r = gem.begin(); r != gem.end(); ++r) {
        const GEMDetId& id = (*r).first;
        const GEMEtaPartition* ep = gemg.etaPartition(id);
        for (auto x = (*r).second.first; x != (*r).second.second; ++x) {
          m << id.region() << id.station() << id.ring() << id.chamber() << id.layer() << id.roll() << x->strip() << x->bx();
          const GlobalPoint g = ep ? ep->toGlobal(ep->centreOfStrip(x->strip())) : GlobalPoint();
          m.xyz(ep ? &g : nullptr, ep ? strips(*ep, ep->specificTopology(), x->strip(), x->strip() + 1) : Box());
        }
      }
    }
    inventory(dir, es, tk, topo, cg);
    std::ofstream fed(f("fed.csv"));
    for (const std::string& n : {evt, std::string("Raw data size per FED: the readout board that sends one data block per event to the DAQ."),
                                std::string("fed_id: FED number (see DataFormats/FEDRawData/interface/FEDNumbering.h)"),
                                std::string("subsystem: the detector the FED belongs to; other = not listed in FEDNumbering.h (588-589)"),
                                std::string("bytes: size of the FED's data block in this event including its header and trailer")})
      fed << "# " << n << '\n';
    fed << "fed_id,subsystem,bytes\n";
    const auto& raw = e.get(rawT_);
    for (int id = 0; id <= FEDNumbering::lastFEDId(); ++id)
      if (const auto n = raw.FEDData(id).size()) fed << id << ',' << fedSystem(id) << ',' << n << '\n';
  }

private:
  // every readout element of the detector, read out or not, with its address and size in channels: the modules,
  // layers and rolls of the tracker and muon system, and the cells of the calorimeters (the ECAL barrel is a full
  // 170 x 360 grid and is not listed). Into inventory/; the same for every event
  void inventory(const std::string& dir, const edm::EventSetup& es, const TrackerGeometry& tk, const TrackerTopology& topo, const CaloGeometry& cg) const {
    const std::string d = dir + "/inventory/";
    std::filesystem::create_directories(d);
    const std::string all = "Every readout element of the CMS detector, whether read out in this event or not.";
    const std::string ends = "Global CMS coordinates in cm.";
    const std::string axes = "ux uy uz / vx vy vz: global direction of the module's local x axis (along which row or strip numbers grow) / local y axis (column numbers)";
    {
      Inv p(d + "pixel.csv", {all, "Pixel modules. " + ends, "part: BPix or FPix; layer: BPix layer 1-4 or FPix disk 1-3",
                              "side: FPix 1 = negative z and 2 = positive z; 0 for BPix", "ladder: BPix ladder or FPix blade",
                              "module: BPix module (1-8 along z) or FPix panel (1-2)", "nrows ncols: pixels per module",
                              "cx cy cz: module centre", axes},
            "detid,part,layer,side,ladder,module,nrows,ncols,cx,cy,cz,ux,uy,uz,vx,vy,vz");
      for (const auto* du : tk.detUnits()) {
        const auto* pu = dynamic_cast<const PixelGeomDetUnit*>(du);
        if (!pu) continue;
        const DetId id = du->geographicalId();
        const bool b = id.subdetId() == 1;
        const auto& s = pu->surface();
        p << id.rawId() << (b ? "BPix" : "FPix") << topo.layer(id) << (b ? 0u : topo.pxfSide(id))
          << (b ? topo.pxbLadder(id) : topo.pxfBlade(id)) << (b ? topo.pxbModule(id) : topo.pxfPanel(id))
          << pu->specificTopology().nrows() << pu->specificTopology().ncolumns();
        p.p(s.position()).v(s.toGlobal(LocalVector(1, 0, 0))).v(s.toGlobal(LocalVector(0, 1, 0))).end();
      }
      Inv q(d + "strip.csv", {all, "Strip sensors. " + ends, "part: TIB TID TOB TEC; layer: barrel layer or disk",
                              "side: TID/TEC 1 = negative z and 2 = positive z; 0 in the barrel", "ring: TID 1-3 or TEC 1-7; 0 in the barrel",
                              "stereo: 1 = stereo sensor of a double-sided module", "nstrips: strips in the sensor", "cx cy cz: sensor centre", axes},
            "detid,part,layer,side,ring,stereo,nstrips,cx,cy,cz,ux,uy,uz,vx,vy,vz");
      const char* part[] = {"", "", "", "TIB", "TID", "TOB", "TEC"};
      for (const auto* du : tk.detUnits()) {
        const auto* su = dynamic_cast<const StripGeomDetUnit*>(du);
        if (!su) continue;
        const DetId id = du->geographicalId();
        const int sd = id.subdetId();
        const auto& s = su->surface();
        q << id.rawId() << part[sd] << topo.layer(id) << (sd == 4 || sd == 6 ? topo.side(id) : 0u)
          << (sd == 4 ? topo.tidRing(id) : sd == 6 ? topo.tecRing(id) : 0u) << int(topo.isStereo(id)) << su->specificTopology().nstrips();
        q.p(s.position()).v(s.toGlobal(LocalVector(1, 0, 0))).v(s.toGlobal(LocalVector(0, 1, 0))).end();
      }
    }
    {
      Inv e(d + "ecal_endcap.csv", {all, "ECAL endcap crystals: zside -1/+1, ix iy 1-100"}, "zside,ix,iy");
      for (const DetId& i : cg.getValidDetIds(DetId::Ecal, EcalEndcap)) {
        const EEDetId id(i);
        e << id.zside() << id.ix();
        e.o << id.iy() << '\n';
      }
      Inv s(d + "ecal_preshower.csv", {all, "ECAL preshower sensors (32 strips each). " + ends,
                                       "x1 y1 / x32 y32: centre of strip 1 / strip 32"},
            "zside,plane,six,siy,x1,y1,x32,y32");
      std::map<std::array<int, 4>, std::array<float, 4>> sensors;
      for (const DetId& i : cg.getValidDetIds(DetId::Ecal, EcalPreshower)) {
        const ESDetId id(i);
        if (id.strip() != 1 && id.strip() != 32) continue;
        const auto g = cg.getGeometry(id)->getPosition();
        auto& v = sensors[{{id.zside(), id.plane(), id.six(), id.siy()}}];
        v[id.strip() == 1 ? 0 : 2] = g.x();
        v[id.strip() == 1 ? 1 : 3] = g.y();
      }
      for (const auto& [k, v] : sensors) {
        s << k[0] << k[1] << k[2] << k[3] << v[0] << v[1] << v[2];
        s.o << v[3] << '\n';
      }
      Inv h(d + "hcal.csv", {all, "HCAL cells (calibration channels excluded): part HB HE HO HF"}, "part,ieta,iphi,depth");
      for (const auto& [sub, name] : std::vector<std::pair<HcalSubdetector, const char*>>{
               {HcalBarrel, "HB"}, {HcalEndcap, "HE"}, {HcalOuter, "HO"}, {HcalForward, "HF"}})
        for (const DetId& i : cg.getValidDetIds(DetId::Hcal, sub)) {
          const HcalDetId id(i);
          h << name << id.ieta() << id.iphi();
          h.o << id.depth() << '\n';
        }
    }
    {
      const auto& dtg = es.getData(dtGT_);
      Inv t(d + "dt.csv", {all, "DT layers. " + ends, "first nwires: first wire number and number of wires",
                           "x0 y0 z0 / x1 y1 z1: midpoint of the first / last wire"},
            "wheel,station,sector,superlayer,layer,first,nwires,x0,y0,z0,x1,y1,z1");
      for (const auto* l : dtg.layers()) {
        const DTLayerId id = l->id();
        const auto& tp = l->specificTopology();
        t << id.wheel() << id.station() << id.sector() << id.superlayer() << id.layer() << tp.firstChannel() << tp.channels();
        t.p(l->toGlobal(LocalPoint(tp.wirePosition(tp.firstChannel()), 0, 0)))
            .p(l->toGlobal(LocalPoint(tp.wirePosition(tp.lastChannel()), 0, 0))).end();
      }
      const auto& cscg = es.getData(cscGT_);
      Inv c(d + "csc.csv", {all, "CSC layers (ME1/1: ring 1 = ME1/1b; ring 4 = ME1/1a). " + ends, "nstrips nwiregroups: per layer",
                            "s1 / sN: centre of strip 1 / the last strip; w1 / wN: centre of wire group 1 / the last one"},
            "zside,station,ring,chamber,layer,nstrips,nwiregroups,s1x,s1y,s1z,sNx,sNy,sNz,w1x,w1y,w1z,wNx,wNy,wNz");
      for (const auto* l : cscg.layers()) {
        const CSCDetId id = l->id();
        const auto* lg = l->geometry();
        const int ns = lg->numberOfStrips(), nw = lg->numberOfWireGroups();
        const auto strip = [&](int k) { return l->toGlobal(lg->topology()->localPosition(MeasurementPoint(k - 0.5, 0))); };
        c << (id.endcap() == 1 ? 1 : -1) << id.station() << id.ring() << id.chamber() << id.layer() << ns << nw;
        c.p(strip(1)).p(strip(ns));
        if (nw > 0) c.p(l->toGlobal(lg->localCenterOfWireGroup(1))).p(l->toGlobal(lg->localCenterOfWireGroup(nw)));
        else c << "" << "" << "" << "" << "" << "";
        c.end();
      }
      const auto& rpcg = es.getData(rpcGT_);
      Inv r(d + "rpc.csv", {all, "RPC rolls. " + ends, "x1 y1 z1 / xN yN zN: centre of strip 1 / the last strip"},
            "region,ring,station,sector,layer,subsector,roll,nstrips,x1,y1,z1,xN,yN,zN");
      for (const auto* roll : rpcg.rolls()) {
        const RPCDetId id = roll->id();
        r << id.region() << id.ring() << id.station() << id.sector() << id.layer() << id.subsector() << id.roll() << roll->nstrips();
        r.p(roll->toGlobal(roll->centreOfStrip(1))).p(roll->toGlobal(roll->centreOfStrip(roll->nstrips()))).end();
      }
      const auto& gemg = es.getData(gemGT_);
      Inv g(d + "gem.csv", {all, "GEM eta partitions. " + ends, "x0 y0 z0 / xN yN zN: centre of strip 0 / the last strip"},
            "region,station,ring,chamber,layer,roll,nstrips,x0,y0,z0,xN,yN,zN");
      for (const auto* ep : gemg.etaPartitions()) {
        const GEMDetId id = ep->id();
        g << id.region() << id.station() << id.ring() << id.chamber() << id.layer() << id.roll() << ep->nstrips();
        g.p(ep->toGlobal(ep->centreOfStrip(0))).p(ep->toGlobal(ep->centreOfStrip(ep->nstrips() - 1))).end();
      }
    }
  }

  const std::string dir_;
  const edm::EDGetTokenT<FEDRawDataCollection> rawT_;
  const edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> pixT_;
  const edm::EDGetTokenT<edm::DetSetVector<SiStripDigi>> stripT_;
  const edm::EDGetTokenT<EBDigiCollection> ebT_;
  const edm::EDGetTokenT<EEDigiCollection> eeT_;
  const edm::EDGetTokenT<ESDigiCollection> esT_;
  const edm::EDGetTokenT<QIE11DigiCollection> hbheT_;
  const edm::EDGetTokenT<QIE10DigiCollection> hfT_;
  const edm::EDGetTokenT<HODigiCollection> hoT_;
  const edm::EDGetTokenT<DTDigiCollection> dtT_;
  const edm::EDGetTokenT<CSCStripDigiCollection> cscsT_;
  const edm::EDGetTokenT<CSCWireDigiCollection> cscwT_;
  const edm::EDGetTokenT<CSCComparatorDigiCollection> csccT_;
  const edm::EDGetTokenT<RPCDigiCollection> rpcT_;
  const edm::EDGetTokenT<GEMDigiCollection> gemT_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkT_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoT_;
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloT_;
  const edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGT_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGT_;
  const edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGT_;
  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGT_;
};

DEFINE_FWK_MODULE(EventCsvDumper);
