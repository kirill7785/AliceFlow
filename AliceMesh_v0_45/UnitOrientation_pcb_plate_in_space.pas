unit UnitOrientation_pcb_plate_in_space;
// pcb плата всегла лежит в плоскости XOY но модель может быть нарисована
// как с номалью к верхней крышке OZ так и с нормалью OY.
// В случае нормали OY предусмотрено преобразование платы XOY->XOZ.
// 29.11.2022

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls;

type
  TFormOrientation_pcb_plate = class(TForm)
    RadioGroupOrientation_pcb: TRadioGroup;
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormOrientation_pcb_plate: TFormOrientation_pcb_plate;

implementation

{$R *.dfm}

end.
