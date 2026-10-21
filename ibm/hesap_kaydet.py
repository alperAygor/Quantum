# IBM Quantum hesabini bu bilgisayara bir kez kaydeder. Anahtar ekranda gorunmez, dosyalara yazilmaz (Qiskit'in kendi ayar dosyasina gider).
from getpass import getpass
from qiskit_ibm_runtime import QiskitRuntimeService
tok = getpass("IBM Quantum API anahtari (yapistir, Enter): ").strip()
crn = input("Instance CRN (IBM Quantum ana sayfasinda 'Instances' altinda; bos birakabilirsin): ").strip() or None
QiskitRuntimeService.save_account(channel="ibm_quantum_platform", token=tok, instance=crn, set_as_default=True, overwrite=True)
svc = QiskitRuntimeService()
print("Kaydedildi. Kullanilabilir donanimlar:", [b.name for b in svc.backends(simulator=False, operational=True)])
