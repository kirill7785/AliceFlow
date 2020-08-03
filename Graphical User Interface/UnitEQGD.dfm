object EGDForm: TEGDForm
  Left = 379
  Top = 129
  AutoSize = True
  Caption = 'Basic parameters'
  ClientHeight = 449
  ClientWidth = 241
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object PanelGlobal: TPanel
    Left = 0
    Top = 0
    Width = 241
    Height = 449
    Color = clMoneyGreen
    TabOrder = 0
    object GBCFlZ: TGroupBox
      Left = 0
      Top = 0
      Width = 217
      Height = 449
      Caption = 'Variables solved'
      TabOrder = 0
      object CBFlow: TCheckBox
        Left = 16
        Top = 23
        Width = 137
        Height = 17
        Caption = 'Flow (velocity/pressure)'
        TabOrder = 0
        OnClick = CBFlowClick
      end
      object RadioGroupFlowRegime: TRadioGroup
        Left = 13
        Top = 126
        Width = 193
        Height = 41
        Caption = 'Flow regime'
        Columns = 2
        ItemIndex = 0
        Items.Strings = (
          'Laminar'
          'Turbulent')
        TabOrder = 1
        OnClick = RadioGroupFlowRegimeClick
      end
      object GroupBoxTurbulentModel: TGroupBox
        Left = 13
        Top = 187
        Width = 193
        Height = 89
        Caption = 'Turbulent model'
        TabOrder = 2
        object ComboBoxturbulentmodel: TComboBox
          Left = 21
          Top = 29
          Width = 169
          Height = 21
          ItemIndex = 0
          TabOrder = 0
          Text = 'Zero Equation Model (RANS)'
          OnChange = ComboBoxturbulentmodelChange
          Items.Strings = (
            'Zero Equation Model (RANS)'
            'Smagorinsky Model (LES)'
            'RNG (LES)'
            'Spalart - Allmares (RANS) [1992]'
            'K-Omega SST (RANS) [1993]')
        end
        object BEditTurb: TButton
          Left = 104
          Top = 56
          Width = 75
          Height = 25
          Caption = 'Edit'
          TabOrder = 1
          Visible = False
          OnClick = BEditTurbClick
        end
      end
      object ButtonidFlow: TButton
        Left = 139
        Top = 421
        Width = 75
        Height = 25
        Caption = 'Apply'
        TabOrder = 3
        OnClick = ButtonidFlowClick
      end
      object GBGravity: TGroupBox
        Left = 13
        Top = 290
        Width = 161
        Height = 125
        Caption = 'gravity vector'
        Color = clMoneyGreen
        ParentColor = False
        TabOrder = 4
        object Lgx: TLabel
          Left = 16
          Top = 32
          Width = 11
          Height = 13
          Caption = 'gx'
        end
        object Lgy: TLabel
          Left = 16
          Top = 64
          Width = 11
          Height = 13
          Caption = 'gy'
        end
        object Lgz: TLabel
          Left = 16
          Top = 96
          Width = 11
          Height = 13
          Caption = 'gz'
        end
        object lblgx: TLabel
          Left = 112
          Top = 24
          Width = 24
          Height = 13
          Caption = 'm/s2'
        end
        object lblgy: TLabel
          Left = 112
          Top = 56
          Width = 24
          Height = 13
          Caption = 'm/s2'
        end
        object lbldz: TLabel
          Left = 112
          Top = 88
          Width = 24
          Height = 13
          Caption = 'm/s2'
        end
        object Egx: TEdit
          Left = 40
          Top = 24
          Width = 65
          Height = 21
          TabOrder = 0
        end
        object Egy: TEdit
          Left = 40
          Top = 56
          Width = 65
          Height = 21
          TabOrder = 1
        end
        object Egz: TEdit
          Left = 40
          Top = 88
          Width = 65
          Height = 21
          TabOrder = 2
        end
      end
      object CheckBoxStaticStructural: TCheckBox
        Left = 21
        Top = 103
        Width = 97
        Height = 17
        Caption = 'Static Structural'
        TabOrder = 5
      end
      object GroupBoxTemperature: TGroupBox
        Left = 13
        Top = 46
        Width = 136
        Height = 51
        Caption = 'Temperature'
        TabOrder = 6
        object ComboBoxTemperature: TComboBox
          Left = 13
          Top = 20
          Width = 108
          Height = 21
          ItemIndex = 1
          TabOrder = 0
          Text = 'Finite Volume Method'
          OnChange = ComboBoxTemperatureChange
          Items.Strings = (
            'none'
            'Finite Volume Method'
            'Finite Element Method'
            'NetWork_T solver')
        end
      end
    end
  end
  object ApplicationEvents1: TApplicationEvents
    OnMessage = ApplicationEvents1Message
    Left = 160
    Top = 31
  end
end
