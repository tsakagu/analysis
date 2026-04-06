#include "TriggerTile.h"

#include <algorithm>  // for sort
#include <numeric>

//______________________________________
void TriggerTile::Reset()
{
  m_tileNumber = 0;
  m_tileEta.clear();
  m_tilePhi.clear();
  m_tileEnergyAdc.clear();
}

//______________________________________
void TriggerTile::identify(std::ostream& out) const
{
  out << " I am a Fired Trigger Tile List" << std::endl;
  out << " There are " << m_tileNumber << " fired trigger tiles" << std::endl;

  for (unsigned int i = 0; i < m_tileNumber; i++)
  {
    std::cout << "Tile " << i << ": "
              << "(eta, phi) = " << "(" << m_tileEta[i] << ", " << m_tilePhi[i] << "), "
              << "E (ADC) = " << m_tileEnergyAdc[i] << std::endl;
  }
}

//______________________________________
void TriggerTile::AddTile(unsigned short ieta, unsigned short iphi, unsigned short eadc)
{
  m_tileNumber++;
  m_tileEta.push_back(ieta);
  m_tilePhi.push_back(iphi);
  m_tileEnergyAdc.push_back(eadc);
}

//______________________________________
void TriggerTile::SortTiles()
{
  std::vector<int> sort_indices(m_tileNumber);
  std::iota(sort_indices.begin(), sort_indices.end(), 0);
  std::sort(sort_indices.begin(), sort_indices.end(), [&](size_t i, size_t j)
            { return m_tileEnergyAdc[i] > m_tileEnergyAdc[j]; });

  std::vector<int> tile_eta_corrected(m_tileNumber);
  std::vector<int> tile_phi_corrected(m_tileNumber);
  std::vector<int> tile_energy_corrected(m_tileNumber);
  for (int i = 0; i < m_tileNumber; i++)
  {
    tile_eta_corrected[i] = m_tileEta[sort_indices[i]];
    tile_phi_corrected[i] = m_tilePhi[sort_indices[i]];
    tile_energy_corrected[i] = m_tileEnergyAdc[sort_indices[i]];
  }
  for (int i = 0; i < m_tileNumber; i++)
  {
    m_tileEta[i] = tile_eta_corrected[i];
    m_tilePhi[i] = tile_phi_corrected[i];
    m_tileEnergyAdc[i] = tile_energy_corrected[i];
  }
}
