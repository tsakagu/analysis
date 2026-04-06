#ifndef TRIGGER_TRIGGERTILE_H
#define TRIGGER_TRIGGERTILE_H

#include <phool/PHObject.h>

#include <vector>

///
class TriggerTile : public PHObject
{
 public:
  TriggerTile() = default;
  virtual ~TriggerTile() override = default;

  /// Clear Event from memory
  void Reset() override;
  void identify(std::ostream& out = std::cout) const override;

  void AddTile(unsigned short seta, unsigned short sphi, unsigned short eadc);

  unsigned short get_tile_eta(int i) const
  {
    if (i < m_tileNumber && i >= 0)
    {
      return m_tileEta[i];
    }
    else
    {
      std::cout << "TriggerTile:: Warning, index out of bounds." << std::endl;
      return -1;
    }
  }

  unsigned short get_tile_phi(int i) const
  {
    if (i < m_tileNumber && i >= 0)
    {
      return m_tilePhi[i];
    }
    else
    {
      std::cout << "TriggerTile:: Warning, index out of bounds." << std::endl;
      return -1;
    }
  }

  unsigned short get_tile_energy_adc(int i) const
  {
    if (i < m_tileNumber && i >= 0)
    {
      return m_tileEnergyAdc[i];
    }
    else
    {
      std::cout << "TriggerTile:: Warning, index out of bounds." << std::endl;
      return -1;
    }
  }

  unsigned short get_tile_number() const
  {
    return m_tileNumber;
  }

  void SortTiles();

 private:  // so the ClassDef does not show up with doc++
  unsigned short m_tileNumber = 0;
  std::vector<unsigned short> m_tileEta;
  std::vector<unsigned short> m_tilePhi;
  std::vector<unsigned short> m_tileEnergyAdc;

  ClassDefOverride(TriggerTile, 1);
};

#endif
