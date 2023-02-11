object Form_debug_panel: TForm_debug_panel
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'Form_debug_panel'
  ClientHeight = 282
  ClientWidth = 546
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox_debug: TGroupBox
    Left = 0
    Top = 0
    Width = 546
    Height = 282
    Caption = 'debug_panel'
    TabOrder = 0
    object Label1: TLabel
      Left = 16
      Top = 27
      Width = 58
      Height = 13
      Caption = 'variable1 = '
    end
    object Edit_free_debug_param: TEdit
      Left = 80
      Top = 24
      Width = 121
      Height = 21
      TabOrder = 0
      Text = '0.0'
    end
  end
end
