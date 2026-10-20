// Progressive enhancement only: native image access and navigation work without JS.
const filter = document.querySelector('#family');
if (filter) {
  filter.closest('.filters').hidden = false;
  filter.addEventListener('change', () => {
    let count = 0;
    document.querySelectorAll('#material-grid .card').forEach(card => {
      card.hidden = filter.value !== 'all' && card.dataset.family !== filter.value;
      if (!card.hidden) count++;
    });
    document.querySelector('#result-count').textContent = `${count} material${count === 1 ? '' : 's'}`;
  });
}
