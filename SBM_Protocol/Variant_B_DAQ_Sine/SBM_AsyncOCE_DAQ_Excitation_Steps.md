# SBM con excitación AsyncOCE (seno generado por el DAQ en ao3)

*24-10-2026. Modo por defecto de `SBM_Protocol`: `Cfg.ExcMode = "sine"`.*

## Qué cambia respecto a las versiones anteriores

| | `gen_puff` / `gen_ao3` (anteriores) | **`sine` (esta versión)** |
|---|---|---|
| Quién genera el seno | Rigol DG922 Pro (ráfaga de 1 ciclo) | **el DAQ, en ao3** |
| Qué sale por ao3 | trigger TTL hacia el Rigol | **el seno, directo al amplificador** |
| Air puff controller | ON (gen_puff) / OFF | **OFF** |
| Filas AO del waveform | 3 (gen_puff) / 4 (gen_ao3) | **4**: X, Y, B-trigger, excitación |
| Sincronismo | re-trigger cada 57 A-lines | **por construcción**: el AO se reloja con el A-line clock (`/Dev1/PFI0`), 57 muestras por período a 285 µs |

**Señal en ao3:** `(ExcOffset_V + ExcVpp_V/2 · sin(2π n/P)) · env(n)`.
- Por defecto **0 a 1 V** (1 Vpp, offset 0.5 V), el mismo nivel que entregaba el Rigol en QH (1 Vpp, 500 mV).
- `env` es una rampa Hann de 5 períodos (`ExcRampPeriods`) al inicio y al final, así que ao3 **arranca y termina en 0 V**, sin escalón sobre el piezo.
- El seno está encendido 50 períodos (≈ 14 ms, `PreamblePeriods`) antes del frame 1 y hasta después del último registro.
- El AsyncOCE original usaba 0 a 0.95 V; si el amplificador espera señal bipolar, poner `ExcOffset_V = 0`.

## Por qué hay que tocar LabVIEW (lo que pasó el 18/19-10)

- Los sidecars de esas adquisiciones muestran que el script node **sí entregó 4 filas** (`size(waveform) = 4 × 81 225`, `Info.LabVIEW.AO_channels = 4`), y aun así ao3 quedó en silencio y DAQmx no dio error.
- Eso significa que `LabVIEW_Setup_Mirrors_and_Btrig.vi` usó la **tarea de 3 canales** (`Dev1/ao0:2`) y descartó la fila 4.
- La tarea de 4 canales (`Dev1/ao0:3` / `Dev1/ao4:5, Dev1/ao2:3`) está dentro de una estructura case con etiqueta **`"AsyncOCE"`**. Con Scan Type = SBM no se entra a ese caso.
- La etiqueta es de tipo string, así que el selector puede venir del *Scan Type* o del *Transducer* del Air puff controller; el archivo no permite distinguirlo. La solución de abajo cubre ambos casos.

---

## Paso 1. MATLAB (PC del Twin)

1. Copiar la carpeta `SBM_Protocol` actualizada (ya la tienes) y verificar en `SBM_TwinOCT_Config.m`:
   ```matlab
   Cfg.ExcMode     = "sine";
   Cfg.ExcVpp_V    = 1.0;    % Vpp en la entrada del amplificador
   Cfg.ExcOffset_V = 0.5;    % 0..1 V (0 si el amplificador es bipolar)
   ```
2. Revisar los vectores antes de LabVIEW:
   ```matlab
   chk = SBM_TwinOCT_PlotVectors;   % debe decir Channels: 4 y "All pass/fail checks PASSED"
   ```
   En la figura, la fila 4 (ao3) debe ser un seno de 285 µs entre 0 y 1 V, con rampa al inicio y al final.
3. Si cambias el FOV (por ejemplo 18 × 18 mm), lo pone LabVIEW (Scan Length → volts con HVC/VVC). Mira `MaxVel_Vpms` y `MaxAcc_Vpms2` en `chk`: 18 mm pide ~1.5 veces la velocidad de galvo de 12 mm.

## Paso 2. LabVIEW: habilitar la tarea de 4 canales para SBM (una vez)

1. Abre `LabVIEW_Setup_Mirrors_and_Btrig.vi` y **guarda una copia de respaldo** (*File → Save As → Copy*).
2. Pulsa **Ctrl+E** (diagrama). Luego **Ctrl+F** → *Text* → escribe `Dev1/ao0:3` → *Find*. LabVIEW salta a la constante string del canal de 4 salidas.
3. Esa constante está dentro de un case interno `"G1", Default` (G1 = canales ao0:3; Default = ao4:5, ao2:3). Sube un nivel: la estructura que lo contiene tiene el selector con **`"AsyncOCE"`**.
4. **Mira qué entra a su terminal selector (?)**: haz clic en el cable y sigue el resaltado hasta la fuente.
   - Si viene de **Scan Type**, basta el paso 5.
   - Si viene del **Air puff controller** (Unbundle → `Enum`/Transducer), haz el paso 5 y además pon Transducer = **SBM** en el preset (Paso 3). El typedef `AirpuffController.ctl` ya tiene el ítem `SBM`.
5. Clic derecho en el selector → **Edit Case Label…** → cambia `"AsyncOCE"` por **`"AsyncOCE", "SBM"`** → OK.
6. Repite la búsqueda **Ctrl+F → `AsyncOCE`**. En este VI hay **otra** estructura case con etiqueta `"AsyncOCE"` (contiene un case `True`). Edítala igual: `"AsyncOCE", "SBM"`. Si encuentras más cases `"AsyncOCE"` en el mismo VI, edítalos todos.
7. Comprueba que la flecha Run no esté rota y guarda (**Ctrl+S**).

> [!note] Por qué es seguro
> El driver AsyncOCE ya entrega 4 filas y pasa por ese mismo caso; ahora SBM usa exactamente el mismo camino. Para los demás modos (OCE, QH…) no cambia nada.

## Paso 3. LabVIEW: preset `SBM_37A_21B` en el master

1. Abre `ATS9360_Master_SSOCT_Trial_3.vi` y pulsa **Ctrl+E**. En el panel frontal, clic derecho sobre el borde del cluster **Air puff controller** → **Find → Local Variables**. Recorre los resultados (doble clic o Ctrl+G) hasta el que está dentro del case **`"SBM_37A_21B"`** del evento `"Set AcqMode Params": Value Change`.
2. En su **Bundle By Name**:

   | Terminal | Constante |
   |---|---|
   | Enable Puff | **False** |
   | Transducer | **SBM**, o **AsyncOCE** si en el Paso 2.4 el selector es el Scan Type (cualquiera sirve en ese caso) |
   | Puff Frequency / ON Length / Delay | da igual (puff apagado); si los agregaste para `gen_puff`, déjalos o bórralos |

3. Revisa también, en el mismo case, que el preset deje Scan Type = **SBM**, A-lines per B-frame **1197**, No. Sync Samples **20**, B-frames **1**, No. Cluster **57**, B-frame Signal Delay **30** y Scan Length según el experimento (12 o 18 mm).
4. Si en el case `"SBM"` del diagrama hay escrituras al *Air puff controller* en un frame **Enabled** de un *Diagram Disable* (lo viste el 20-10), elimínalas o haz que coincidan con esta tabla.
5. Guarda. Al usarlo, selecciona otro preset y luego `SBM_37A_21B` otra vez para que el evento se dispare.

## Paso 4. Cableado

1. **Desconecta el Rigol** del amplificador, o deja su salida en OFF.
2. Conecta el **BNC de ao3 → entrada del amplificador del piezo** (donde antes entraba la salida CH1 del Rigol). El BNC de ao3 ya no va a la entrada de trigger del Rigol.
3. El AO del DAQ entrega ±10 V con pocos mA. Una entrada de amplificador de alta impedancia está bien; **no la cargues con 50 Ω**.

## Paso 5. Verificación con osciloscopio (antes de adquirir)

Pon el CH1 del osciloscopio en ao3 (en T con la entrada del amplificador). Trigger en flanco de subida a 0.5 V, Single.

| Escala | Debe verse |
|---|---|
| 100 µs/div | seno de **285.0 µs** (*Measure → Period*) entre **0 y 1 V**; escalonado fino (57 muestras/período, a 200 kHz) |
| 50 ms/div | un tren de seno de ~**0.41 s** con rampa de entrada y salida (~1.4 ms), y 0 V antes y después |

- **ao3 plano en 0 V** → LabVIEW sigue usando 3 canales: revisa el Paso 2 (etiqueta del case y selector) y el Paso 3 (Transducer / Enable Puff).
- **Error DAQmx** "number of channels in the data does not match the task" → filas ≠ canales: comprueba que `Cfg.ExcMode = "sine"` (4 filas) y que entró el caso de 4 canales.
- Opcional: `SBM_TwinOCT_CheckExcitation` con el CSV del osciloscopio comprueba la frecuencia, el lock y la amplitud.

## Paso 6. Adquirir y comprobar

1. Adquiere una prueba de 285 µs (nombre según el log).
2. En MATLAB (en el PC del Twin o en el de procesamiento):
   ```matlab
   SBM_CheckAcquisitionConsistency(binFile)   % debe dar ALL CONSISTENT y "AO rows ... 4 4"
   ```
3. `Main_Processing_SingleExp_SBM_EMF` debe imprimir **`EXCITATION DETECTED`**. El sidecar debe decir `ExcMode = "sine"`.

## Volver a QH (MB)

Selecciona el preset `QH`; su case reescribe el puff (400 Hz / 0.5 ms / 0.1 ms). Reconecta ao3 → trigger del Rigol y la salida del Rigol → amplificador. La edición `"AsyncOCE", "SBM"` no afecta a QH.
